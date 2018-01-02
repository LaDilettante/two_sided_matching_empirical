rm(list = ls())

source("match2sided.R")
library("tidyverse")
library("coda")
time_start <- Sys.time()

# ---- Data ----
dat <- read.table("../data_clean/gss18cat.raw", header=T)

n_i <- nrow(dat) # number_of_employees, Job acceptances (elements in 1:n_j)
n_j <- length(unique(dat$occ17))    # number of employers, includes unemployment
p_i <- 4 # number of worker characteristics per worker; including the intercept
# in this example it measures quality of employee
p_j <- 2  # number of job characteristics per job; in this example, job quality

# ---- Prepare choice ----
choice <- dat$occ17 + 1

# ---- Prepare opp, the opportunity set ----
opp <- matrix(F, n_i, n_j)  # The opportunity matrix T=offer,F=no offer
opp[cbind(1:n_i, choice)] <- T  # people are offered jobs they have!
opp[, 1] <- T                     # Unemployment always offered

# ---- Prepare ww ----
ww <- matrix(NA, n_j, p_j)
for (i in 1:n_j) {
  ww[i, 1] <- unique(dat$presmean[choice == i])
  ww[i, 2] <- unique(dat$autmean[choice == i])
}
ww <- scale(ww) # rescale (for better numerics )
colnames(ww) <- c("prestige", "autonomy")
rownames(ww) <- c("Unemployment", "Professionals, Self employed",
                  "Professionals, Salaried",
                  "Managers", "Salesmen, Other", "Proprietors",
                  "Clerical", "Salesmen, Retail",
                  "Craftsmen, Manufacturing",
                  "Craftsmen, Other", "Craftsmen, Construction",
                  "Operatives, Manufacturing", "Operatives, Other",
                  "Service", "Laborers, Manufacturing",
                  "Laborers, Other", "Farmers", "Farm laborers")

# ---- Prepare xx ----
# matrix nworkers x nx of worker characteristics
# including column of ones for an intercept
one <- rep(1, n_i)
xx <- scale(dat[, c("educ", "age", "nonwhite")])
xx <- cbind(one, xx)

# ---- Run MCMC ----

res <- match2sided(iter = 20000,
                  C_alpha = (0.4 ** 2) * diag(ncol(ww)), 
                  C_beta = (0.025 ** 2) * diag(ncol(xx)), 
                  frac_opp = 0.25,
                  ww = ww, xx = xx,
                  choice = choice, opp = opp)
print(res$acceptance_rate)
plot(res$lp[, 1], type='l',
     xlab = 'iteration', ylab = 'log joint pdf')
plot(res$lp[, 2], type='l', xlab = 'iteration', ylab = 'lp_A')
plot(res$lp[, 3], type='l', xlab = 'iteration', ylab = 'lp_O')
saveRDS(res, file = paste0("../result/labor-", 
                           format(Sys.time(), "%Y-%m-%d-%H%M"), ".RData"))

# ---- Results and Diagnostics ----

WARMUP <- res$mcmc_settings$iter / 5

alpha <- mcmc(res$alpha) %>% window(start = WARMUP)
plot(alpha)
summary(alpha)

beta_educ <- mcmc(res$beta[, 'educ', ])
plot(beta_educ[, c('Professionals, Salaried', 'Farm laborers')])
# plot(beta_educ)

betastar_educ <- mcmc(res$betastar[, 2, ])
plot(betastar_educ[, c(3, 18)])

beta_age <- mcmc(res$beta[, 'age', ])
plot(beta_age[, c('Professionals, Salaried', 'Farm laborers')])
# plot(beta_age)

plot(density(res$beta[WARMUP:dim(res$beta)[1], 1, 1]))
for (i in 2:dim(res$beta)[2]) {
  lines(density(res$beta[WARMUP:dim(res$beta)[1], i, 1]))
}

plot(density(res$beta[WARMUP:dim(res$beta)[1], 2, 1]))
for (i in 3:dim(res$beta)[2]) {
  lines(density(res$beta[WARMUP:dim(res$beta)[1], i, 1]))
}

boxplot(educ ~ occ17, data = dat %>% filter(occ17 %in% c(2, 3, 10, 17)),
        names = c("professional", "managers",
                  "construction worker", "farm worker"),
        ylab = "years of education",
        main = "sample statistics of education")

par(mfrow = c(1, 2))
boxplot(educ ~ occ17, data = dat %>% filter(occ17 %in% c(2, 17)),
        names = c("professional", "farm worker"),
        ylab = "years of education",
        main = "sample statistics of education")

plot(density(beta_educ[ , c("Professionals, Salaried")]),
     xlim = c(-0.5, 0.5),
     main = "estimated firms' preference for education")
lines(density(beta_educ[ , c("Farm laborers")]), lty = 2)
legend(-0.45, 6, c("profressional", "farm worker"), lty = c(1, 2))
par(mfrow = c(1, 1))

# mu_beta
plot(mcmc(res$mu_beta))

# ---- Regress beta ----

X_beta <- dat %>% 
  dplyr::select(occ17, presmean, autmean) %>%
  distinct() %>% arrange(occ17)
X_beta$occ <- c("Unemployment", "Professionals, Self employed",
                "Professionals, Salaried",
                "Managers", "Salesmen, Other", "Proprietors",
                "Clerical", "Salesmen, Retail",
                "Craftsmen, Manufacturing",
                "Craftsmen, Other", "Craftsmen, Construction",
                "Operatives, Manufacturing", "Operatives, Other",
                "Service", "Laborers, Manufacturing",
                "Laborers, Other", "Farmers", "Farm laborers")

posterior_mean_beta <- data.frame(t(apply(res$beta, MARGIN = c(2, 3), mean))) %>%
  mutate(occ = rownames(.))
dat_beta <- inner_join(X_beta, posterior_mean_beta, by = "occ")

summary(lm(educ ~ presmean + autmean, data = dat_beta))
summary(lm(age ~ presmean + autmean, data = dat_beta))
summary(lm(nonwhite ~ presmean + autmean, data = dat_beta))

# ---- Data viz ---

viz_dat <- dat
viz_dat <- viz_dat %>%
  mutate(college = educ > 12)

ggplot(data = viz_dat) +
  geom_density(aes(presmean)) +
  facet_wrap(~college, labeller = label_both) +
  labs(x = "prestige of the job")
