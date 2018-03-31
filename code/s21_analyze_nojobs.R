rm(list = ls())
library("tidyverse")
library("reticulate")
library("coda")
source("match2sided.R")
source("0_plot_functions.R")

use_python("/Users/anh/miniconda3/bin/python")

set.seed(1)
py_run_file("s11_labor_nojobs.py")
py$choice <- py$choice + 1
py$xx[, 2] <- py$xx[, 2] - mean(py$xx[, 2])
py$xx[, 3] <- py$xx[, 3] - mean(py$xx[, 3])
#py$ww <- scale(py$ww, center = TRUE, scale = FALSE)

# ---- Estimating beta from true_opp ----

# one-sided beta
sapply(c(2, 3, 4, 5, 6), function(j) {
  m <- glm(py$true_opp[, j] ~ . - 1, family = binomial(link = "logit"),
      data = as.data.frame(py$xx))
  coef(m)
})

# true beta
t(py$true_beta)

# ---- Estimating alpha from true_opp ----

# Only look at people who got the opportunity set c(1, 0, 1, 1)
idx <- which(colSums(t(py$true_opp) == c(1, 0, 0, 1, 1, 1)) == 6)
choice_tmp <- py$choice[idx]
choice_tmp[py$choice[idx] == 4] <- 2
choice_tmp[py$choice[idx] == 5] <- 3
choice_tmp[py$choice[idx] == 6] <- 4

# Calculate by hand
cl_nllik <- function(alpha) {
  ww_tmp <- py$ww[c(1, 4, 5, 6), ]
  wa <- ww_tmp %*% alpha
  lse_wa <- log(sum(exp(wa)))
  - sum(wa[choice_tmp] - lse_wa)
}
# one-sided alpha (not very good, cuz of high correlation)
optim(c(0, 0), cl_nllik)

# true alpha
py$true_alpha

# ---- Check the effect of initial starting_opp values ----

# starting_opp <- matrix(NA, nrow = dim(jobs)[1], ncol = dim(jobs)[2])
# starting_opp[, 1] <- 1
# for (i in 2:ncol(starting_opp)) {
#   probs <- mean(true_opp[, i])
#   starting_opp[, i] <- sample(c(0, 1), size = nrow(starting_opp),
#                               prob = c(1 - probs, probs), replace = TRUE)
# }
# starting_opp[cbind(1:nrow(starting_opp), choice)] <- 1
# 
# summary(glm(starting_opp[, 2] ~ . - 1, family = binomial(link = "logit"),
#             data = as.data.frame(xx)))
# summary(glm(starting_opp[, 3] ~ . - 1, family = binomial(link = "logit"),
#             data = as.data.frame(xx)))
# summary(glm(starting_opp[, 4] ~ . - 1, family = binomial(link = "logit"),
#             data = as.data.frame(xx)))


# ---- MCMC ----

iter <- 1e4
thin <- 5

prior <- list(alpha = list(mu = rep(0, py$p_j), Tau = solve(diag(rep(100, py$p_j)))),
              mu_beta = list(mu = rep(0, py$p_i),
                             Tau = solve(diag(rep(100, py$p_i)))),
              Tau_beta = list(nu = py$p_i + 2,
                              Sinv = solve(diag(rep(100, py$p_i))))) # E(\Sigma) = Sinv (Hoff p110)

start_time <- Sys.time()
starting_alpha <- rnorm(py$p_j, mean = 0, sd = 1)
starting_beta <- matrix(rnorm(p_i * n_j, 0, sd = 1), p_i, n_j)
res <- match2sided(iter = iter, t0 = iter + 1,
                   C_alpha = diag(c(0.01, 0.1) ** 2), 
                   C_beta = diag(c(0.01, rep(0.005, py$p_i - 1)) ** 2),
                   starting_alpha = rep(0, py$p_j),
                   starting_beta = matrix(0, nrow = py$p_i, ncol = py$n_j),
                   frac_opp = 0.25, prior = prior,
                   ww = py$ww, xx = py$xx,
                   choice = py$choice, starting_opp = py$obs_opp,
                   reserve_choice = TRUE,
                   to_save = c("alpha", "beta", "opp"),
                   file = paste("../result/sim_nojobs_", start_time), write = FALSE)
cat("Sim labor no jobs done\n")
cat("Running time:", difftime(Sys.time(), start_time, units = "mins"), "mins")
cat("Memory used: ", gc()[2, 2])
colMeans(res$ok)

start <- iter / 2

pdf(paste0("../figure/alpha ", start_time, ".pdf"), w = 7.5, h = 5)
par(oma=c(0, 0, 3, 0))
my.plot.mcmc(mcmc(res$alpha[, ]), parameters = py$true_alpha)
mtext("Workers' preference param", side = 3, line = 0, outer = TRUE)
dev.off()

mcmcse::mcse.multi(res$beta[start:iter, , 2]) # beta for the 1st employer
pdf(paste0("../figure/beta_emp1 ", start_time, ".pdf"), 2 = 7.5, h = 7.5)
par(oma=c(0, 0, 3, 0))
my.plot.mcmc(mcmc(res$beta[, , 2]), parameters = c(-3.613, py$true_beta[2, 2:4]))
mtext("Employer 1's preference params", side = 3, line = 0, outer = TRUE)
dev.off()

my.plot.mcmc(mcmc(res$beta[, , 2]), parameters = c(-9.0, 0.2, 0.2))
my.plot.mcmc(mcmc(res$betastar[, , 2]), parameters = py$true_beta[2, ])

mcmcse::mcse.multi(res$beta[start:iter, , 3]) # beta for the 2nd employer
pdf(paste0("../figure/beta_emp2 ", start_time, ".pdf"), 2 = 7.5, h = 7.5)
par(oma=c(0, 0, 3, 0))
my.plot.mcmc(mcmc(res$beta[, , 3]), parameters = py$true_beta[3, ])
mtext("Employer 2's preference params", side = 3, line = 0, outer = TRUE)
dev.off()

mcmcse::mcse.multi(res$beta[start:iter, , 4]) # beta for the 3rd employer
par(oma=c(0, 0, 3, 0))
my.plot.mcmc(mcmc(res$beta[, , 4]), parameters = py$true_beta[4, ])
mtext("Employer 3's preference params", side = 3, line = 0, outer = TRUE)

my.plot.mcmc(mcmc(res$beta[, , 5]), parameters = py$true_beta[5, ])
my.plot.mcmc(mcmc(res$beta[, , 6]), parameters = py$true_beta[6, ])

beta_2 <- res$beta[, 2, ] # the beta for the second employee's covariate (educ)
my.plot.mcmc(mcmc(beta_2))
colMeans(beta_2[5000:10000, ])

# Check updating of opp
for (i in 1:5) {
  heatmap2(res$opp[i, , ])
}

# See the percentage of offfer for each employers and compare with true_opp
# This could be an approximation for how stringent the
tmp <- t(apply(res$opp, 1, 
               function(sampled_opp) colMeans(sampled_opp) - colMeans(true_opp)))
par(mfrow = c(4, 2))
par(oma=c(0, 0, 3, 0))
f_plot_mcmc(tmp[seq(1, iter, thin), 2], main = "Employer 2")
f_plot_mcmc(tmp[seq(1, iter, thin), 3], main = "Employer 3")
f_plot_mcmc(tmp[seq(1, iter, thin), 4], main = "Employer 4")
f_plot_mcmc(tmp[seq(1, iter, thin), 5], main = "Employer 5")
mtext("Difference in Offer Rate b/w Sampled Opp Set and True Opp Set", 
      side = 3, line = 0, outer = TRUE)

# Average difference between the sampled_opp and true_opp for each employer
tmp <- t(apply(res$opp, 1, 
               function(sampled_opp) colMeans(abs(sampled_opp - true_opp))))
par(mfrow = c(4, 2))
par(oma=c(0, 0, 3, 0))
f_plot_mcmc(tmp[seq(1, iter, thin), 2], main = "Employer 2")
f_plot_mcmc(tmp[seq(1, iter, thin), 3], main = "Employer 3")
f_plot_mcmc(tmp[seq(1, iter, thin), 4], main = "Employer 4")
f_plot_mcmc(tmp[seq(1, iter, thin), 5], main = "Employer 5")
mtext("Avg Difference b/w Sampled Opp Set and True Opp Set", 
      side = 3, line = 0, outer = TRUE)

# ---- Correlation between opp and beta ----

offer_rate <- t(apply(res$opp, 1, colMeans))
intercept <- res$beta[ , 1, ] # intercept

par(mfrow = c(2, 2))
plot(offer_rate[, 2], intercept[, 2])
plot(offer_rate[, 3], intercept[, 3])
plot(offer_rate[, 4], intercept[, 4])

cor(offer_rate, intercept)