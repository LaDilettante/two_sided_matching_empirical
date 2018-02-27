rm(list = ls())
library(RcppCNPy)
library(coda)
library(zoo) # moving average for MCMC
source("match2sided2.R")
source("0_functions.R")

ww <- npyLoad("ww.npy", "numeric")
xx <- npyLoad("xx.npy", "numeric")
choice <- npyLoad("choice.npy") + 1 # Python is index 0, R is 1
jobs <- npyLoad("jobs.npy")
dim(jobs) <- c(dim(jobs), 1) # Turn jobs into a 3D array
job_utilities <- npyLoad("job_utilities.npy")
true_opp <- npyLoad("true_opp.npy")
obs_opp <- npyLoad("obs_opp.npy")

# ---- Estimating beta from true_opp ----

summary(glm(true_opp[, 2] ~ . - 1, family = binomial(link = "logit"),
        data = as.data.frame(xx)))
summary(glm(true_opp[, 3] ~ . - 1, family = binomial(link = "logit"),
            data = as.data.frame(xx)))
summary(glm(true_opp[, 4] ~ . - 1, family = binomial(link = "logit"),
            data = as.data.frame(xx)))

# ---- Estimating alpha from true_opp ----

# Only look at people who got the opportunity set c(1, 0, 1, 1)
idx <- which(colSums(t(true_opp) == c(1, 1, 0, 1)) == 4)
jobs2 <- jobs[idx, c(T, T, F, T), , drop = FALSE]
choice2 <- choice[idx]
choice2[choice[idx] == 4] <- 3

# Calculate by hand
cl_nllik <- function(alpha) {
  job_utils <- apply(jobs2, c(1, 2), function(x) alpha %*% x)
  wa <- job_utils[cbind(1:nrow(job_utils), choice2)]
  lse_wa <- apply(job_utils, 1, function(x) log(sum(exp(x))))
  - sum(wa - lse_wa)
}
optimize(cl_nllik, lower = -10, upper = 10)

# ---- Check the effect of initial starting_opp values ----

starting_opp <- matrix(NA, nrow = dim(jobs)[1], ncol = dim(jobs)[2])
starting_opp[, 1] <- 1
for (i in 2:ncol(starting_opp)) {
  probs <- mean(true_opp[, i])
  starting_opp[, i] <- sample(c(0, 1), size = nrow(starting_opp),
                              prob = c(1 - probs, probs), replace = TRUE)
}
starting_opp[cbind(1:nrow(starting_opp), choice)] <- 1

summary(glm(starting_opp[, 2] ~ . - 1, family = binomial(link = "logit"),
            data = as.data.frame(xx)))
summary(glm(starting_opp[, 3] ~ . - 1, family = binomial(link = "logit"),
            data = as.data.frame(xx)))
summary(glm(starting_opp[, 4] ~ . - 1, family = binomial(link = "logit"),
            data = as.data.frame(xx)))


# ---- MCMC ----

iter <- 1e4
start <- iter / 2
thin <- 1
p_i <- ncol(xx)
prior <- list(alpha = list(mu = 0, Tau = solve(matrix(0.01))),
              mu_beta = list(mu = rep(0, p_i),
                             Tau = solve(diag(rep(0.01, p_i)))),
              Tau_beta = list(nu = p_i + 2,
                              Sinv = solve(diag(rep(0.01, p_i)))))

start_time <- Sys.time()
res <- match2sided(iter = iter, t0 = iter + 1,
                   C_alpha = matrix(0.01), 
                   C_beta = diag(c(0.01, rep(0.005, ncol(xx) - 1)) ** 2),
                   starting_alpha = c(0),
                   frac_opp = 0.25, prior = prior,
                   jobs = jobs, xx = xx,
                   choice = choice, opp = obs_opp,
                   to_save = c("alpha", "beta", "opp"),
                   file = paste("../result/sim", start_time), write = FALSE)

pdf(paste0("../figure/alpha ", start_time, ".pdf"), w = 7.5, h = 5)
par(oma=c(0, 0, 3, 0))
my.plot.mcmc(mcmc(res$alpha), parameters = c(0.05))
mtext("Workers' preference param", side = 3, line = 0, outer = TRUE)
dev.off()

plot(mcmc(res$alphastar))

colMeans(res$alpha[start:iter, , drop = FALSE])

mcmcse::mcse.multi(res$alpha)

mcmcse::mcse.multi(res$beta[start:iter, , 2]) # beta for the 1st employer
pdf(paste0("../figure/beta_emp1 ", start_time, ".pdf"), 2 = 7.5, h = 7.5)
par(oma=c(0, 0, 3, 0))
my.plot.mcmc(mcmc(res$beta[, , 2]), parameters = c(-9.0, 0.2, 0.2))
mtext("Employer 1's preference params",  side = 3, line = 0, outer = TRUE)
dev.off()

my.plot.mcmc(mcmc(res$beta[, , 2]), parameters = c(-9.0, 0.2, 0.2))
plot(mcmc(res$betastar[, , 2]))

mcmcse::mcse.multi(res$beta[start:iter, , 3]) # beta for the 2nd employer
pdf(paste0("../figure/beta_emp2 ", start_time, ".pdf"), 2 = 7.5, h = 7.5)
par(oma=c(0, 0, 3, 0))
my.plot.mcmc(mcmc(res$beta[, , 3]), parameters = c(-10.0, 0.1, 0.2))
mtext("Employer 2's preference params", side = 3, line = 0, outer = TRUE)
dev.off()

mcmcse::mcse.multi(res$beta[start:iter, , 4]) # beta for the 3rd employer
par(oma=c(0, 0, 3, 0))
my.plot.mcmc(mcmc(res$beta[, , 4]), parameters = c(-5.0, 0.5, -0.05))
mtext("Employer 3's preference params", side = 3, line = 0, outer = TRUE)

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
