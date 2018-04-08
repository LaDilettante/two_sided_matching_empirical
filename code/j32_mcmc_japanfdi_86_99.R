rm(list = ls())

library("tidyverse")
library("coda")
source("match2sided.R")

# ---- Load MNC data ----

df_mnc_raw <- readRDS("../data_clean/japanFDI_86_99.RData")

# Cut down the data
# Only manufacturing, founded in 1996
C_SIC <- "Manufacturing"
C_FDTN_YR <- 1996
df_mnc <- df_mnc_raw %>%
  filter(sic_brd == C_SIC) %>%
  filter(fdtn_yr == C_FDTN_YR, year >= C_FDTN_YR) %>%
  select(sort_key, ltemp, luscptl, int_r_d, int_exp, nation, year) %>% na.omit() %>%
  # Only keep the data closest to the year 1996
  group_by(sort_key) %>% arrange(year) %>% filter(row_number() == 1) %>%
  ungroup() %>%
  # Only keep countries with at least 2 plants
  group_by(nation) %>% filter(n() >= 5) %>% ungroup() %>%
  filter(nation %in% c("China", "Indonesia", "Malaysia", "Philippines", 
                       "Singapore", "South Korea", "Taiwan", "Thailand",
                       "Vietnam")) %>% # Only keep countries in East / SE Asia
  select(ltemp, luscptl, int_r_d, int_exp, nation)

# Assign numeric id to countries
df_nation_id <- df_mnc %>% select(nation) %>% distinct() %>% 
  arrange(nation) %>%
  mutate(nation_id = 1:n(), 
         iso3c = countrycode::countrycode(nation, "country.name", "iso3c"))

# Merge nation_id back, finalize
df_mnc <- df_mnc %>%
  inner_join(df_nation_id, by = 'nation')

xx <- df_mnc %>%
  select(ltemp, luscptl, int_r_d, int_exp) %>%
  mutate_at(vars(ltemp, luscptl), scale, center = TRUE, scale = FALSE) %>%
  mutate_at(vars(int_r_d, int_exp), scale, center = TRUE, scale = TRUE) %>%
  mutate(one = 1) %>%
  select(one, everything()) %>%
  as.matrix()

# ---- Load country data ----
df_country_raw <- readRDS("../data_clean/country_86_99.RData")
df_country <- df_country_raw %>%
  filter(year == C_FDTN_YR) %>%
  inner_join(df_nation_id, by = c('iso3c')) %>%
  select(country, nation, nation_id, lpop, lgdp, lgdppc, 
         hc) %>%
  # Make sure that the nations are arranged according to nation_id, buz we'll turn em into matrix
  arrange(nation_id) %>% select(-nation_id, -country)

ww <- df_country %>%
  select(lgdp, lgdppc, hc) %>%
  as.matrix()
rownames(ww) <- df_country$nation

# ---- Prepare obs_opp ----

n_i <- nrow(xx) ; p_i <- ncol(xx)
n_j <- nrow(ww) ; p_j <- ncol(ww)

choice <- df_mnc$nation_id
obs_opp <- matrix(FALSE, n_i, n_j)  # The opportunity matrix T=offer,F=no offer
obs_opp[cbind(1:n_i, choice)] <- TRUE  # firms are offered the countries they are in!

# ---- Summary statistics ----

df_table <- df_mnc %>% count(nation) %>% 
  mutate(percent = n / sum(n) * 100)
print(xtable::xtable(df_table),
      floating = FALSE,
      include.rownames = FALSE, file = "../table/japan96_num_obs.tex")

# ---- MCMC ----

# iter <- 4e5, thin <- 10 is often what we run eventually
iter <- 4e5
thin <- 10
prior <- list(alpha = list(mu = rep(0, p_j), Tau = solve(diag(rep(100, p_j)))),
              mu_beta = list(mu = c(0, rep(0, p_i - 1)),
                             Tau = solve(diag(c(1, rep(10, p_i - 1))))),
              Tau_beta = list(nu = 25, # p_i + const
                              Sinv = solve(diag(c(5, rep(10, p_i - 1))))))

start_time <- Sys.time()

starting_alpha <- rep(0, p_j)
starting_beta <- matrix(0, nrow = p_i, ncol = n_j)
res <- match2sided(iter = iter, t0 = iter + 1, thin = thin,
                   C_alpha = diag(c(0.25, 0.25, 0.25) ** 2), 
                   C_beta = diag(c(0.1, 0.1, 0.1, 0.1, 0.1) ** 2),
                   # C_beta_est = C_beta_est,
                   starting_alpha = starting_alpha,
                   starting_beta = starting_beta,
                   frac_opp = 0.15, prior = prior,
                   ww = ww, xx = xx,
                   choice = choice, starting_opp = obs_opp,
                   reserve_choice = FALSE,
                   to_save = c("alpha", "beta", "opp"),
                   file = paste("../result/sim_nojobs_", start_time), write = FALSE)
cat("Japan 96 done\n")
cat("Running time:", difftime(Sys.time(), start_time, units = "mins"), "mins")
cat("Memory used: ", gc()[2, 2])
colMeans(res$ok)

# standardize xx, kinda diffuse prior
paste0("../result/japan96_", start_time, ".RData")
saveRDS(res, paste0("../result/japan96_", start_time, ".RData"))

# ---- Result and Diagnostics ----

pdf(paste0("../figure/japan_post_dens_", start_time, ".pdf"), w = 7, h = 7)
par(mfrow = c(2, 2))
plot(res$lp[, 1], type='l',
     xlab = 'iteration', ylab = 'log posterior density')
plot(res$lp[, 2], type='l', xlab = 'iteration', ylab = 'lp_A')
plot(res$lp[, 3], type='l', xlab = 'iteration', ylab = 'lp_O')
par(mfrow = c(1, 1))
dev.off()



plot(mcmc(res$alpha))
beta_one <- mcmc(res$beta[, 'one', ])
plot(beta_one)

beta_temp <- mcmc(res$beta[, 'ltemp', ])
plot(beta_temp)

plot(mcmc(res$betastar[, 'ltemp', ]))

beta_uscptl <- mcmc(res$beta[, 'luscptl', ])
plot(beta_uscptl)

beta_rd <- mcmc(res$beta[, 'int_r_d', ])
plot(beta_rd)

beta_exp <- mcmc(res$beta[, 'int_exp', ])
plot(beta_exp)


# ---- Investigating opp ----

tmp <- apply(res$opp, c(1), function(mat) {
  colMeans(mat)
})
par(mfrow = c(5, 2))
for (i in 1:n_j) {
  plot(tmp[i, ], type = "l", ylab = i)
}

par(mfrow = c(4, 1))
plot(res$lp[, 1], type='l',
     xlab = 'iteration', ylab = 'log posterior density')
plot(tmp[3, ], type = "l")
plot(res$beta[, 'one', 3], type = "l")
plot(res$beta[, 'ltemp', 3], type = "l")

# Look at sampled firms in Vietnam / or China
id <- df_nation_id %>% filter(nation == "Vietnam") %>% .[["nation_id"]]
mnc_vietnam_obs <- df_mnc %>% 
  filter(nation_id == id) %>%
  colMeans()

# Look at the true offered firms in Vietnam / or China

opp_vietnam <- res$opp[, , id]
mnc_vietnam <- t(apply(opp_vietnam, 1, function(offer) colMeans(xx[which(offer), ])))

plot(mnc_vietnam[, 'ltemp'], type = 'l')
abline(h = mnc_vietnam_obs['ltemp'], col = 'red')

plot(density(xx[, "ltemp"]))
abline(v = mnc_vietnam_obs['ltemp'], col = 'red')
