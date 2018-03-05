rm(list = ls())

source("0_functions.R")
source("match2sided.R")
library("tidyverse")
library("coda")

# ---- Data ----
df_employer_raw <- read_csv("../data_clean/labor_employer_occ5.csv")
ww <- df_employer_raw %>%
  select(pres, aut) %>%
  as.matrix()
rownames(ww) <- df_employer_raw$occ5

df_employee_raw <- read_csv("../data_clean/labor_employee.csv")
xx <- df_employee_raw %>%
  select(educ, age, nonwhite) %>%
  mutate(one = 1) %>%
  select(one, everything()) %>%
  as.matrix()

n_i <- nrow(xx) # number_of_employees, Job acceptances (elements in 1:n_j)
n_j <- nrow(ww)    # number of employers, includes unemployment
p_i <- ncol(xx) # number of worker characteristics per worker; including the intercept
p_j <- ncol(ww)  # number of job characteristics per job; in this example, job quality

# ---- Prepare choice ----
choice <- df_employee_raw$occ5_num + 1

# ---- Prepare opp, the opportunity set ----
obs_opp <- matrix(FALSE, n_i, n_j)  # The opportunity matrix T=offer,F=no offer
obs_opp[cbind(1:n_i, choice)] <- TRUE  # people are offered jobs they have!
obs_opp[, 1] <- TRUE                     # Unemployment always offered

# ---- Run MCMC ----

iter <- 2e5
thin <- 10
prior <- list(alpha = list(mu = rep(0, p_j), Tau = solve(diag(rep(0.01, p_j)))),
              mu_beta = list(mu = rep(0, p_i),
                             Tau = solve(diag(rep(0.01, p_i)))),
              Tau_beta = list(nu = p_i + 2,
                              Sinv = solve(diag(rep(0.01, p_i)))))

start_time <- Sys.time()
res <- match2sided(iter = iter, t0 = iter + 1, thin = thin,
                   C_alpha = diag(c(0.001, 0.005)), 
                   C_beta = diag(c(0.05, 0.005, 0.0025, 0.01) ** 2),
                   starting_alpha = rep(0, p_j),
                   frac_opp = 0.25, prior = prior,
                   ww = ww, xx = xx,
                   choice = choice, opp = obs_opp,
                   to_save = c("alpha", "beta"),
                   file = paste("../result/sim_nojobs_", start_time), write = FALSE)
cat("labor_occ5 done \n")

saveRDS(res, paste0("../result/labor_occ5_", start_time, ".pdf"))

# ---- MCMC analysis ----

colMeans(res$ok)

# ----- Results -----

pdf(paste0("../figure/labor_occ5_post_dens_", start_time, ".pdf"), w = 7, h = 7)
par(mfrow = c(2, 2))
plot(res$lp[, 1], type='l',
     xlab = 'iteration', ylab = 'log posterior density')
plot(res$lp[, 2], type='l', xlab = 'iteration', ylab = 'lp_A')
plot(res$lp[, 3], type='l', xlab = 'iteration', ylab = 'lp_O')
par(mfrow = c(1, 1))
dev.off()

# alpha
pdf(paste0("../figure/labor_occ5_alpha_", start_time, ".pdf"), w = 7, h = 7)
par(oma=c(0, 0, 3, 0))
plot(mcmc(res$alpha))
mtext("Workers' preference parameters", side = 3, line = 0, outer = TRUE)
dev.off()

# beta
beta_educ <- mcmc(res$beta[, 'educ', ])
pdf(paste0("../figure/labor_occ5_beta_educ_", start_time, ".pdf"), w = 7, h = 7)
plot(beta_educ[ , 2:ncol(beta_educ)])
dev.off()

beta_age <- mcmc(res$beta[, 'age', ])
pdf(paste0("../figure/labor_occ5_beta_age_", start_time, ".pdf"), w = 7, h = 7)
plot(beta_age[, 2:ncol(beta_age)])
dev.off()

tmp <- res$beta
dim(tmp) <- c(dim(tmp)[1], dim(tmp)[2] * dim(tmp)[3])
pdf("../figure/correlation_across_beta_adaptive.pdf", w = 7, h = 7.4)
lattice::levelplot(cor(tmp),
                   main = "Correlation across 72 beta estimates\n18 employers, 4 beta each")
dev.off()