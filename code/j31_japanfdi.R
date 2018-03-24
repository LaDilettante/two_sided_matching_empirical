rm(list = ls())

library("tidyverse")
library("coda")
source("match2sided.R")

# ---- Data ----
df_country_raw <- read_csv("../data_clean/japan_country.csv")
ww <- df_country_raw %>%
  select(lgdp, lgdppc, democracy) %>%
  as.matrix()
rownames(ww) <- df_country_raw$country

df_mnc_raw <- read_csv("../data_clean/japan_mnc.csv",
                       col_types = cols(lemp = "d", lcap = "d", 
                                        country = "c", country_id = "i"))
xx <- df_mnc_raw %>%
  select(lemp, lcap) %>%
  mutate(one = 1) %>%
  select(one, everything()) %>%
  as.matrix()

n_i <- nrow(xx) # number_of_employees, Job acceptances (elements in 1:n_j)
n_j <- nrow(ww)    # number of employers, includes unemployment
p_i <- ncol(xx) # number of worker characteristics per worker; including the intercept
p_j <- ncol(ww)  # number of job characteristics per job; in this example, job quality

# ---- Prepare choice ----
choice <- df_mnc_raw$country_id

# ---- Prepare opp ----
# Opportunity sets
obs_opp <- matrix(F, n_i, n_j)  # The opportunity matrix T=offer,F=no offer
obs_opp[cbind(1:n_i, choice)] <- TRUE  # firms are offered the countries they are in!

# ---- Run MCMC ----

iter <- 5e4
thin <- 10
prior <- list(alpha = list(mu = rep(0, p_j), Tau = solve(diag(rep(0.01, p_j)))),
              mu_beta = list(mu = rep(0, p_i),
                             Tau = solve(diag(rep(0.01, p_i)))),
              Tau_beta = list(nu = p_i + 2,
                              Sinv = solve(diag(rep(0.01, p_i)))))

start_time <- Sys.time()
res <- match2sided(iter = iter, t0 = iter + 1, thin = thin,
                   C_alpha = diag(c(0.001, 0.005, 0.0025)), 
                   C_beta = diag(c(0.001, 0.001, 0.001) ** 2),
                   starting_alpha = rep(0, p_j),
                   frac_opp = 0.1, prior = prior,
                   ww = ww, xx = xx,
                   choice = choice, opp = obs_opp,
                   to_save = c("alpha", "beta"),
                   file = paste("../result/sim_nojobs_", start_time), write = FALSE)
cat("japan FDI done \n")
colMeans(res$ok)
  
saveRDS(res, file = paste0("../result/japan-", start_time, ".RData"))

# ---- Results and Diagnostics ----

pdf(paste0("../figure/japan_post_dens_", start_time, ".pdf"), w = 7, h = 7)
par(mfrow = c(2, 2))
plot(res$lp[, 1], type='l',
     xlab = 'iteration', ylab = 'log posterior density')
plot(res$lp[, 2], type='l', xlab = 'iteration', ylab = 'lp_A')
plot(res$lp[, 3], type='l', xlab = 'iteration', ylab = 'lp_O')
par(mfrow = c(1, 1))
dev.off()

# alpha
pdf(paste0("../figure/japan_alpha_", start_time, ".pdf"), w = 7, h = 7)
par(oma=c(0, 0, 3, 0))
plot(mcmc(res$alpha))
mtext("MNCs' preference parameters", side = 3, line = 0, outer = TRUE)
dev.off()

# beta
beta_intercept <- mcmc(res$beta[, 'one', ])
plot(beta_intercept)

beta_temp <- mcmc(res$beta[, 'lemp', ])
pdf(paste0("../figure/japan_beta_educ_", start_time, ".pdf"), w = 7, h = 7)
plot(beta_temp)
dev.off()

beta_lcap <- mcmc(res$beta[, 'lcap', ])
pdf(paste0("../figure/japan_beta_age_", start_time, ".pdf"), w = 7, h = 7)
plot(beta_lcap)
dev.off()

tmp <- res$beta
dim(tmp) <- c(dim(tmp)[1], dim(tmp)[2] * dim(tmp)[3])
pdf("../figure/japan_correlation_across_beta_adaptive.pdf", w = 7, h = 7.4)
lattice::levelplot(cor(tmp),
                   main = "Correlation across 72 beta estimates\n18 employers, 4 beta each")
dev.off()