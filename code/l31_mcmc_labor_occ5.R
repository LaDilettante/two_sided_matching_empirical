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
  # Scale to mean 0
  mutate_at(vars(educ, age), scale, center = TRUE, scale = FALSE) %>% 
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

iter <- 2e4
thin <- 10
prior <- list(alpha = list(mu = rep(0, p_j), Tau = solve(diag(rep(100, p_j)))),
              mu_beta = list(mu = rep(0, p_i),
                             Tau = solve(diag(rep(100, p_i)))),
              Tau_beta = list(nu = p_i + 2,
                              Sinv = solve(diag(rep(100, p_i)))))

start_time <- Sys.time()
res <- match2sided(iter = iter, t0 = iter + 1, thin = thin,
                   C_alpha = diag(c(0.01, 0.05) ** 2), 
                   C_beta = diag(c(0.05, 0.01, 0.01, 0.05) ** 2),
                   starting_alpha = rep(0, p_j),
                   starting_beta = matrix(0, nrow = p_i, ncol = n_j),
                   frac_opp = 0.25, prior = prior,
                   ww = ww, xx = xx,
                   choice = choice, starting_opp = obs_opp,
                   reserve_choice = TRUE,
                   to_save = c("alpha", "beta", "opp"),
                   file = paste("../result/sim_nojobs_", start_time), write = FALSE)
cat("labor_occ5 done \n")
cat("Running time:", difftime(Sys.time(), start_time, units = "mins"), "mins")
cat("Memory used: ", gc()[2, 2])
colMeans(res$ok)


# standardize xx, 1st
paste0("../result/labor_occ5_", start_time, ".RData")
saveRDS(res, paste0("../result/labor_occ5_", start_time, ".RData"))

# ----- Results -----

par(mfrow = c(2, 2))
plot(res$lp[, 1], type='l',
     xlab = 'iteration', ylab = 'log posterior density')
plot(res$lp[, 2], type='l', xlab = 'iteration', ylab = 'lp_A')
plot(res$lp[, 3], type='l', xlab = 'iteration', ylab = 'lp_O')
par(mfrow = c(1, 1))

# alpha
plot(mcmc(res$alpha))

# beta
beta_one <- mcmc(res$beta[, 'one', ])
plot(beta_one[, 2:ncol(beta_one)])

beta_educ <- mcmc(res$beta[, 'educ', ])
plot(beta_educ[ , 2:ncol(beta_educ)])

beta_age <- mcmc(res$beta[, 'age', ])
plot(beta_age[, 2:ncol(beta_age)])

beta_nonwhite <- mcmc(res$beta[, 'nonwhite', ])
plot(beta_nonwhite)

tmp <- res$beta
dim(tmp) <- c(dim(tmp)[1], dim(tmp)[2] * dim(tmp)[3])
pdf("../figure/correlation_across_beta_adaptive.pdf", w = 7, h = 7.4)
lattice::levelplot(cor(tmp),
                   main = "Correlation across 72 beta estimates\n18 employers, 4 beta each")
dev.off()