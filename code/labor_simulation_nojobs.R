rm(list = ls())
library(RcppCNPy)
library(coda)
library(zoo) # moving average for MCMC
source("match2sided.R")
source("0_functions.R")

ww <- npyLoad("ww_nojobs.npy", "numeric")
xx <- npyLoad("xx_nojobs.npy", "numeric")
choice <- npyLoad("choice_nojobs.npy") + 1 # Python is index 0, R is 1
true_opp <- npyLoad("true_opp_nojobs.npy", dotranspose = FALSE)
obs_opp <- npyLoad("obs_opp_nojobs.npy")

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
ww2 <- ww[c(T, T, F, T)]
choice2 <- choice[idx]
choice2[choice[idx] == 4] <- 3

# Calculate by hand
cl_nllik <- function(alpha) {
  wa <- ww2 * alpha
  lse_wa <- log(sum(exp(wa)))
  - sum(wa[choice2] - lse_wa)
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

iter <- 1e5
p_i <- ncol(xx)
p_j <- ncol(ww)
prior <- list(alpha = list(mu = rep(0, p_j), Tau = solve(diag(rep(0.01, p_j)))),
              mu_beta = list(mu = rep(0, p_i),
                             Tau = solve(diag(rep(0.01, p_i)))),
              Tau_beta = list(nu = p_i + 2,
                              Sinv = solve(diag(rep(0.01, p_i)))))

res <- match2sided(iter = iter, t0 = iter + 1,
                   C_alpha = diag(c(0.01, 0.1) ** 2), 
                   C_beta = diag(c(0.1, rep(0.01, ncol(xx) - 1)) ** 2),
                   frac_opp = 0.25, prior = prior,
                   ww = ww, xx = xx,
                   choice = choice, opp = obs_opp)

plot(mcmc(res$alpha))
