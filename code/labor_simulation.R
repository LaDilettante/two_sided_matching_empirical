rm(list = ls())
library(RcppCNPy)
library(coda)
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

iter <- 10000
p_i <- ncol(xx)
prior <- list(alpha = list(mu = 0, Tau = matrix(0.01)),
              mu_beta = list(mu = rnorm(p_i),
                             Tau = solve(diag(rep(0.01, p_i)))),
              Tau_beta = list(nu = p_i + 2,
                              Sinv = solve(diag(rep(0.01, p_i)))))



res <- match2sided(iter = iter, t0 = iter + 1,
                   C_alpha = matrix(0.01), 
                   C_beta = (0.01 ** 2) * diag(ncol(xx)),
                   starting_alpha = c(0),
                   frac_opp = 0.25, prior = prior,
                   jobs = jobs, xx = xx,
                   choice = choice, opp = starting_opp)

start <- iter / 2
plot(mcmc(res$alpha))
plot(mcmc(res$alphastar))
colMeans(res$alpha[start:iter, , drop = FALSE])

mcmcse::mcse.multi(res$alpha)

mcse.multi(res$beta[start:iter, , 2]) # beta for the 2nd employer
plot(mcmc(res$beta[, , 2]))
plot(mcmc(res$betastar[, , 2]))

mcse.multi(res$beta[start:iter, , 4]) # beta for the 3rd employer
plot(mcmc(res$beta[, , 3]))

beta_2 <- res$beta[, 2, ] # the beta for the second employee's covariate (educ)
plot(mcmc(beta_2))
colMeans(beta_2[5000:10000, ])

# Check updating of opp
for (i in 1:5) {
  heatmap2(res$opp[i, , ])
}

tmp <- t(apply(res$opp, 1, colMeans))
tmp <- tmp - colMeans(true_opp)
plot(tmp[, 2], type = "l")
plot(tmp[, 3], type = "l")
plot(tmp[, 4], type = "l")

sum(abs(res$opp[2, , ] - true_opp))
