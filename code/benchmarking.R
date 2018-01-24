library("profvis")
library("microbenchmark")

# ---- Fixing replicate being slow ----

n_i <- 2000
n_j <- 20
frac_opp <- 0.25

sample(n_i * n_j, size = n_i * floor(frac_opp * n_j), replace = FALSE)

microbenchmark(
  replicate(n_i, sample(2:n_j, size=floor(frac_opp * n_j), replace = FALSE)),
  sample(n_i * n_j, size = n_i * floor(frac_opp * n_j), replace = FALSE)
)

microbenchmark(
  replicate(n_i, sample(2:n_j, size=floor(frac_opp * n_j), replace = FALSE))
)


# ---- Checking the entire process ----

res <- match2sided(iter = 20000,
                   C_alpha = (0.4 ** 2) * diag(ncol(ww)), 
                   C_beta = (0.025 ** 2) * diag(ncol(xx)),
                   frac_opp = 0.25,
                   ww = ww, xx = xx,
                   choice = choice, opp = opp)

calculate_C <- function(X_sample, t, C, Xbar, eps = 0.01) {
  d <- length(X_sample)
  sd <- (2.4 ** 2) / d
  
  C_old <- C
  Xbar_old <- Xbar
  
  Xbar <- 1 / (t - 1) * ((t - 2) * Xbar_old + X_sample)
  C <- (t - 3) / (t - 2) * C_old +
    sd / (t - 2) * ((t - 2) * Xbar_old %*% t(Xbar_old) -
                      (t - 1) * Xbar %*% t(Xbar) +
                      X_sample %*% t(X_sample) +
                      eps * diag(d))
  
  return(list(C = C, Xbar = Xbar))
}

iter <- 1000
t0 <- 100
C_alpha = (0.4 ** 2) * diag(ncol(ww))
C_beta = (0.025 ** 2) * diag(ncol(xx))
frac_opp = 0.25
profvis({
  n_i <- nrow(xx)
  n_j <- nrow(ww)
  p_i <- ncol(xx)
  p_j <- ncol(ww)
  eps <- 0.01 # Small constant for beta proposal
  
  prior <- list(alpha = list(mu = rnorm(p_j),
                             Tau = solve(diag(abs(rnorm(p_j))))),
                mu_beta = list(mu = rnorm(p_i),
                               Tau = solve(diag(abs(rnorm(p_i))))),
                Tau_beta = list(nu = p_i + 2,
                                Sinv = solve(diag(abs(rnorm(p_i))))))
  
  # ---- Starting values ----
  alpha <- rep(0, p_j)            # worker preferences  
  exp_WA <- exp(ww %*% alpha) # linear predictors of
  pA_den <- opp %*% exp_WA # vector of denominators in p(A | O, alpha)
  
  # beta starting values (from 1-sided logit estimates)

  beta <- matrix(0, p_i, n_j)   # employer preferences
  for(j in 1:n_j) {
    y <- as.numeric(opp[, j])
    mod <- glm(y ~ . - 1, family=binomial,
               data=as.data.frame(xx) )
    beta[, j] <- mod$coef
  }  
  betastar <- beta
  
  XB <- xx %*% beta # worker side linear predictors (big matrix), n_i x n_j, same as opp
  
  # mu, Sigma starting values
  mu_beta <- prior$mu_beta$mu
  Tau_beta <- prior$Tau_beta$Sinv
  
  # ---- Pre compute ----
  tmp <- as.matrix(ww[choice, ])
  wa <- apply(tmp, 2, sum) # sum of characteristics of accepted jobs; used in alpha update
  
  # ---- Initialize storage ----
  acceptance_rate <- rep(0, 3)                  # Metropolis acceptance rates
  asave <- matrix(NA, iter, p_j)   # saved alphas
  bsave <- array(NA, dim = c(iter, p_i, n_j)) # saved betas
  bstarsave <- array(NA, dim = c(iter, p_i, n_j)) # saved betas proposal
  mu_betasave <- matrix(NA, iter, p_i)
  Tau_betasave <- array(NA, dim = c(iter, p_i, p_i))
  logpost <- matrix(NA, iter, 3) # posterior density, i.e. joint, lp_A, lp_O
  timesave <- matrix(NA, iter, 3)
  
  # Store results
  logpost[1, 1] <- joint_lpdf(opp, choice,
                              alpha, beta, mu_beta, Tau_beta, prior,
                              ww, xx)
  logpost[1, 2] <- f_logp_A(opp, choice, alpha, ww)
  logpost[1, 3] <- f_logp_O(opp, beta, xx)
  asave[1, ] <- alpha
  bsave[1, , ] <- beta # vectorize
  bstarsave[1, , ] <- betastar
  mu_betasave[1, ] <- mu_beta
  Tau_betasave[1, ,] <- Tau_beta
  timesave[1, ] <- c(0, 0, 0)
  
  # ---- Loop ----
  start100 <- Sys.time()
  for (i in 2:iter) {
    start <- Sys.time()
    # Update opp
    start_opp <- Sys.time()
    
    new <- sample(n_i * n_j, size = n_i* floor(frac_opp * n_j), replace = FALSE)
    accepted_job <- 1:n_i + (choice - 1) * n_i # note that matrix index is by column as default
    new <- setdiff(new, accepted_job) # remove accepted job from flipping consideration
    unemployment <- 1:n_i # unemp is the first column
    new <- setdiff(new, unemployment) # remove unemployment from flipping consideration
    ind <- arrayInd(new, .dim = c(n_i, n_j))
    
    my_logmh_opp <- logmh_opp(opp, new, alpha, beta, ww, xx)
    ok_opp <- log(runif(n_i)) <= my_logmh_opp
    
    if (any(ok_opp)) {
      old_opp <- opp
      opp[ind] <- !(old_opp[ind]) # Update the opportunity set
      opp[!ok_opp, ] <- old_opp[!ok_opp, ] # Un-update the rows that should not be updated 
      acceptance_rate[1] <- acceptance_rate[1] + mean(ok_opp)
    }
    timesave[i, 1] <- Sys.time() - start_opp
    
    # Update alpha and beta
    start_ab <- Sys.time()
    sd <- (2.4 ** 2) / (p_j + p_i * n_j)
    if (i <= t0) {
      C_ab_est <- as.matrix(Matrix::bdiag(c(list(C_alpha), 
                                            rep(list(C_beta), n_j))))
    } else if (i > t0) {
      if (i == t0 + 1) {
        # Calculate Cs and Xbars for the first time
        alpha_samples <- asave[1:(i - 1), ]
        beta_samples <- bsave[1:(i - 1), , ]
        dim(beta_samples) <- c(i - 1, p_i * n_j)
        ab_samples <- cbind(alpha_samples, beta_samples)
        C_ab_est <- sd * cov(ab_samples) + sd * eps * diag(p_j + p_i * n_j)
        Xbar_ab_est <- colMeans(ab_samples)
      } else {
        # Calculate Cs and Xbars recursively
        # (notice we're passing in old values of Cs and Xbars)
        
        res <- calculate_C(X_sample = c(asave[i - 1, ], bsave[i - 1, , ]), 
                           t = i, 
                           C = C_ab_est, Xbar = Xbar_ab_est)
        C_ab_est <- res$C
        Xbar_ab_est <- res$Xbar
      }
    }
    
    deviation <- c(rmvnorm(1, sigma = C_ab_est))
    alphastar <- alpha + deviation[1:p_j]
    beta_deviation <- matrix(deviation[(p_j + 1):(p_j + p_i * n_j)], 
                             nrow = p_i, ncol = n_j)
    betastar <- beta + beta_deviation
    
    my_logmh_alpha <- logmh_alpha(alpha, alphastar, ww, opp, wa, prior)
    my_logmh_beta <- logmh_beta(beta, betastar, xx, opp, mu_beta, Tau_beta)
    ok_ab <- ifelse(log(runif(1)) <= my_logmh_alpha + my_logmh_beta, T, F)
    if (ok_ab) {
      alpha <- alphastar
      beta <- betastar
      acceptance_rate[2] <- acceptance_rate[2] + 1
      acceptance_rate[3] <- acceptance_rate[3] + 1
    }
    timesave[i, 2] <- Sys.time() - start_ab
    
    # Update mu (multivariate normal)
    mu_beta_posterior <- cond_mu_beta(Tau_beta, prior, beta)
    mu_beta <- mvtnorm::rmvnorm(1, mu_beta_posterior$m,
                                mu_beta_posterior$V)
    mu_beta <- c(mu_beta) # Flatten into a vector
    
    # Update Tau (Wishart)
    Tau_beta_posterior <- cond_Tau_beta(mu_beta, prior, beta)
    Tau_beta <- MCMCpack::rwish(Tau_beta_posterior$nu,
                                Tau_beta_posterior$Sinv)
    
    # Store results
    logpost[i, 1] <- joint_lpdf(opp, choice,
                                alpha, beta, mu_beta, Tau_beta, prior,
                                ww, xx)
    logpost[i, 2] <- f_logp_A(opp, choice, alpha, ww)
    logpost[i, 3] <- f_logp_O(opp, beta, xx)
    asave[i, ] <- alpha
    bsave[i, , ] <- beta # vectorize
    bstarsave[i, , ] <- betastar
    mu_betasave[i, ] <- mu_beta
    Tau_betasave[i, ,] <- Tau_beta
    
    timesave[i, 3] <- Sys.time() - start
    # Increment
    if (i %% 100 == 0) {
      cat("Iteration", i, "done", Sys.time() - start100, "\n")
      start100 <- Sys.time()
    }
  }
  
  if (!is.null(colnames(ww))) colnames(asave) <- colnames(ww)
  if (!is.null(colnames(xx))) {
    dimnames(bsave)[[2]] <- colnames(xx)
    dimnames(bstarsave)[[2]] <- colnames(xx)
    colnames(mu_betasave) <- colnames(xx)
    dimnames(Tau_betasave)[[2]] <- colnames(xx)
    dimnames(Tau_betasave)[[3]] <- rownames(xx)
  }
  if (!is.null(rownames(ww))) {
    dimnames(bsave)[[3]] <- rownames(ww)
    dimnames(bstarsave)[[3]] <- rownames(ww)
  }  
})


microbenchmark(
  replicate(n_i, sample(2:n_j, size=floor(frac_opp * n_j), replace = FALSE)),
  bsave[1:1000, , ],
  bsave[1:10000, , ]
)

microbenchmark(
  alpha_samples <- asave[1:(i - 1), ],
  beta_samples <- bsave[1:(i - 1), , ],
  dim(beta_samples) <- c(i - 1, p_i * n_j),
  ab_samples <- cbind(alpha_samples, beta_samples)  
)

microbenchmark(
  bsave[1:10000, , ],
  bsave[1:10, , ]
)
