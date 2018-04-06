rm(list = ls())
library("tidyverse")
library("reticulate")
library("coda")
source("match2sided.R")
source("0_plot_functions.R")

use_python("/Users/anh/miniconda3/bin/python")

set.seed(1)
py_run_file("s11_labor_nojobs.py")

idx_nooffer <- which(rowSums(py$true_opp) == 0)
py$true_opp <- py$true_opp[-idx_nooffer, ]
py$obs_opp <- py$obs_opp[-idx_nooffer, ]
py$choice <- py$choice[-idx_nooffer]
py$xx <- py$xx[-idx_nooffer, ]
py$n_i <- nrow(py$xx)
py$xx[, 2] <- py$xx[, 2] - mean(py$xx[, 2])
py$xx[, 3] <- py$xx[, 3] - mean(py$xx[, 3])

# Re-calculate true beta after standardizing xx
py$true_beta

py$true_beta_std <- py$true_beta
py$true_beta_std[, 1] <- py$true_beta[, 1] + py$true_beta[, c(2, 3)] %*% colMeans(py$df_xx[-idx_nooffer, c(2, 3)])

# ---- Estimating beta from true_opp ----

# one-sided beta
sapply(c(1:py$n_j), function(j) {
  m <- glm(py$true_opp[, j] ~ . - 1, family = binomial(link = "logit"),
      data = as.data.frame(py$xx))
  coef(m)
})

# true beta
t(py$true_beta_std)

# ---- Estimating alpha from true_opp ----

# Only look at people who got the opportunity set c(1, 0, 1, 1)
idx <- which(colSums(t(py$true_opp) == c(1, 1, 1, 1, 1)) == 5)
df_count <- data.frame(table(py$choice[idx]),
                       py$df_ww)
glm(Freq ~ pres + aut, data = df_count, family = "poisson")

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
# ---- Run ----
start_time <- Sys.time()
starting_alpha <- rep(0, py$p_j)
starting_beta <- matrix(0, py$p_i, py$n_j)
res <- match2sided(iter = iter, t0 = iter + 1,
                   C_alpha = diag(c(0.015, 0.15) ** 2), 
                   C_beta = diag(c(0.01, 0.005, 0.005, 0.01) ** 2),
                   starting_alpha = starting_alpha,
                   starting_beta = starting_beta,
                   frac_opp = 0.25, prior = prior,
                   ww = py$ww, xx = py$xx,
                   choice = py$choice, starting_opp = py$obs_opp,
                   reserve_choice = FALSE,
                   to_save = c("alpha", "beta", "opp"),
                   file = paste("../result/sim_nojobs_", start_time), write = FALSE)
cat("Sim labor no jobs done\n")
cat("Running time:", difftime(Sys.time(), start_time, units = "mins"), "mins")
cat("Memory used: ", gc()[2, 2])
colMeans(res$ok)

# no unemployment, informative prior
saveRDS(res, paste0("../result/simlabor_nojobs_", start_time, ".RData"))

# ---- Analysis ----

start <- iter / 2

my.plot.mcmc(mcmc(res$alpha[, ]), parameters = py$true_alpha)
my.plot.mcmc(mcmc(res$beta[, , 1]), parameters = py$true_beta_std[1, ])
my.plot.mcmc(mcmc(res$beta[, , 2]), parameters = py$true_beta_std[2, ])
my.plot.mcmc(mcmc(res$beta[, , 3]), parameters = py$true_beta_std[3, ])
my.plot.mcmc(mcmc(res$beta[, , 4]), parameters = py$true_beta_std[4, ])
my.plot.mcmc(mcmc(res$beta[, , 5]), parameters = py$true_beta_std[5, ])

my.plot.mcmc(mcmc(res$beta[, 2, 5] / res$beta[, 3, 5]),
             parameters = py$true_beta_std[5, 2] / py$true_beta_std[5, 3])

j <- 4
tmp <- res$beta[start:iter, 2, j] / res$beta[start:iter, 3, j]
true <- py$true_beta_std[j, 2] / py$true_beta_std[j, 3]
q <- quantile(tmp, c(0.05, 0.95))
tmp2 <- tmp[tmp >= q[1] & tmp <= q[2]]
plot(density(tmp2)) ; abline(v = true, col = 'red')


# Check updating of opp
for (i in 1:5) {
  heatmap2(res$opp[i, , ])
}

# See the percentage of offfer for each employers and compare with true_opp
# This could be an approximation for how stringent the employer is
tmp <- t(apply(res$opp, 1, 
               function(sampled_opp) colMeans(sampled_opp)))
plot(tmp[, 1], type = 'l') ; abline(h = colMeans(py$true_opp)[1], col = 'red')
plot(tmp[, 2], type = 'l') ; abline(h = colMeans(py$true_opp)[2], col = 'red')
plot(tmp[, 3], type = 'l') ; abline(h = colMeans(py$true_opp)[3], col = 'red')

# ---- Correlation between opp and beta ----

offer_rate <- t(apply(res$opp, 1, colMeans))
intercept <- res$beta[ , 1, ] # intercept

par(mfrow = c(2, 1))
plot(offer_rate[, 2], type = 'l') ; abline(h = colMeans(py$true_opp)[2], col = 'red')
plot(intercept[, 2], type = 'l') ; abline(h = py$true_beta_std[2, 1], col = 'red')

par(mfrow = c(2, 1))
plot(offer_rate[, 1], type = 'l') ; abline(h = colMeans(py$true_opp)[1], col = 'red')
plot(offer_rate[, 2], type = 'l') ; abline(h = colMeans(py$true_opp)[2], col = 'red')


par(mfrow = c(2, 2))
plot(offer_rate[, 1], intercept[, 1])
plot(offer_rate[, 2], intercept[, 2])
plot(offer_rate[, 3], intercept[, 3])
plot(offer_rate[, 4], intercept[, 4])

cor(offer_rate, intercept)
