rm(list = ls())
ww2 <- matrix(c(18, 40.9, 33.5, 
                0.256, 0.065, 0.123), nrow = 3) # Covariates, taken from real data
true_alpha <- c(0.01, 1) # Assuming we know the true preference parameter
num_obs <- 1000 # Number of observations
num_choices <- 3 # Number of choices
S <- 1000 # Number of simulations
res <- matrix(NA, nrow = S, ncol = length(true_alpha)) # Storage for simulation result
for (s in 1:S) {
  Uij <- matrix(NA, nrow = num_obs, ncol = num_choices)
  for (j in 1:num_choices) {
    Uij[, j] <- sum(true_alpha * ww2[j, ]) + evd::rgumbel(num_obs)
  }
  y <- max.col(Uij)
  cl_nllik <- function(alpha) {
    wa <- c(ww2 %*% alpha)
    lse_wa <- log(sum(exp(wa)))
    - sum( wa[y] - lse_wa )
  }
  res[s, ] <- optim(c(0, 0), cl_nllik)$par
}
colMeans(res)
apply(res, 2, sd) # Really poor precision!

# ---- 1D preference ----

rm(list = ls())
ww2 <- matrix(c(18, 40.9, 33.5), nrow = 3) # Covariates, taken from real data
true_alpha <- c(0.01) # Assuming we know the true preference parameter
num_obs <- 1000 # Number of observations
num_choices <- 3 # Number of choices
S <- 1000 # Number of simulations
res <- matrix(NA, nrow = S, ncol = length(true_alpha)) # Storage for simulation result
for (s in 1:S) {
  Uij <- matrix(NA, nrow = num_obs, ncol = num_choices)
  for (i in 1:num_choices) {
    Uij[, i] <- sum(true_alpha * ww2[i, ]) + evd::rgumbel(num_obs)
  }
  y <- max.col(Uij)
  cl_nllik <- function(alpha) {
    wa <- c(ww2 %*% alpha)
    lse_wa <- log(sum(exp(wa)))
    - sum( wa[y] - lse_wa )
  }
  res[s, ] <- optimize(cl_nllik, lower = -10, upper = 10)$minimum
}
colMeans(res)
apply(res, 2, sd) # Really poor precision!

# ---- Another 1D preference ----

ww2 <- matrix(c(0.256, 0.065, 0.123), nrow = 3) # Covariates, taken from real data
true_alpha <- c(1) # Assuming we know the true preference parameter
num_obs <- 1000 # Number of observations
num_choices <- 3 # Number of choices
S <- 1000 # Number of simulations
res <- matrix(NA, nrow = S, ncol = length(true_alpha)) # Storage for simulation result
for (s in 1:S) {
  Uij <- matrix(NA, nrow = num_obs, ncol = num_choices)
  for (i in 1:num_choices) {
    Uij[, i] <- sum(true_alpha * ww2[i, ]) + evd::rgumbel(num_obs)
  }
  y <- max.col(Uij)
  cl_nllik <- function(alpha) {
    wa <- c(ww2 %*% alpha)
    lse_wa <- log(sum(exp(wa)))
    - sum( wa[y] - lse_wa )
  }
  res[s, ] <- optimize(cl_nllik, lower = -10, upper = 10)$minimum
}
colMeans(res)
apply(res, 2, sd) # Really poor precision!
