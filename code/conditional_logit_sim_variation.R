rm(list = ls())

# Conditional logit where the covariates vary across both choices AND respondents
ww2 <- matrix(c(18, 40.9, 33.5, 0.256, 0.065, 0.123), nrow = 3) # Covariates, taken from real data
true_alpha <- c(0.07, 1) # Assuming we know the true preference parameter
num_obs <- 1000 # Number of observations
num_choices <- 3 # Number of choices
S <- 1000 # Number of simulations
res <- matrix(NA, nrow = S, ncol = length(true_alpha)) # Storage for simulation result
for (s in 1:S) {
  Uij <- matrix(NA, nrow = num_obs, ncol = num_choices)
  # Add some noise so that the covariates vary by respondent as well
  ww2_list <- vector("list", length = num_obs)
  for (i in 1:num_obs) {
    ww2_list[[i]] <- ww2 + mvtnorm::rmvnorm(3, sigma = diag(c(10, 0.05)))
  }
  # Calculate the utility of respondent i for choice j
  for (i in 1:num_obs) {
    for (j in 1:num_choices) {
      Uij[i, j] <- sum(true_alpha * ww2_list[[i]][j, ]) + evd::rgumbel(1)  
    }
  }
  # Respondents choosing what gives the highest utility
  y <- max.col(Uij)
  # Negative log likelihood to be minimized
  cl_nllik <- function(alpha) {
    ll <- 0
    for (i in 1:num_obs) {
      wa <- c(ww2_list[[i]] %*% alpha)
      lse_wa <- log(sum(exp(wa)))
      ll <- ll + wa[y[i]] - lse_wa
    }
    -ll
  }
  res[s, ] <- optim(c(0, 0), cl_nllik)$par
}
colMeans(res)
apply(res, 2, sd)