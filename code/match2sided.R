library("mvtnorm")
library("plyr")
Rcpp::sourceCpp("0_functions.cpp")
Rcpp::sourceCpp("match2sided.cpp")

#' Two sided matching model
#'
#' @param opp A matrix of opportunity
#' @param choice A vector of observed choice
#' @param alpha Multinomial preference
#' @param beta Binomial preference
#' @param mu_beta Hyperparameter of beta
#' @param sigma_beta Hyperparameter of beta
#' @param prior A list of prior parameters
NULL

# ---- MH and conditionals ----

#' Calculate the log joint pdf, useful for testing
#'
#' @rdname match2sided
#' @importFrom mvtnorm mvtnorm::dmvnorm
f_logp_A <- function(opp, choice, alpha, ww) {
  w <- as.matrix(ww[choice, ])
  logp_A <- sum(w %*% alpha) - sum(log(opp %*% exp(ww %*% alpha)))
  return(logp_A)
}
f_logp_O <- function(opp, beta, xx) {
  XB <- xx %*% beta
  logp_O <- sum(opp * XB) - sum(log1p(exp(XB)))
  return(logp_O)
}

joint_lpdf <- function(opp, choice, alpha, beta,
                       mu_beta, Tau_beta,
                       prior, ww, xx) {
  n_i <- length(choice)
  logp_A <- f_logp_A(opp, choice, alpha, ww)
  # w <- as.matrix(ww[choice, ])
  # logp_A <- sum(w %*% alpha) - sum(log(opp %*% exp(ww %*% alpha)))
  
  logp_O <- f_logp_O(opp, beta, xx)
  # XB <- xx %*% beta
  # logp_O <- sum(opp * XB) - sum(log1p(exp(XB)))
  return(logp_A + logp_O +
           mvtnorm::dmvnorm(alpha, prior$alpha$mu, solve(prior$alpha$Tau), log = TRUE) +
           sum(mvtnorm::dmvnorm(t(beta), mu_beta, solve(Tau_beta), log = TRUE)) +
           mvtnorm::dmvnorm(mu_beta, prior$mu_beta$mu, solve(prior$mu_beta$Tau), log = TRUE) +
           log(MCMCpack::dwish(Tau_beta, prior$Tau_beta$nu, prior$Tau_beta$S)))
}

f_pA_den <- function(opp, ww, alpha) {
  exp_WA <- exp(ww %*% alpha)
  return(c(opp %*% exp_WA))
}

#' log_mh opp_set
#' @param new Indexes of jobs to be flipped, n_i x num_new_offers_each_i
logmh_opp <- function(opp, new, alpha, beta, ww, xx) {
  # browser()
  n_i <- nrow(xx)
  num_new_offers <- length(new) / n_i
  n_j <- nrow(ww)
  
  indnew <- cbind(rep(1:n_i, each = num_new_offers), new)
  oo <- opp[indnew] # current offer status of the jobs under consideration
  plusminus <- ifelse(oo, -1, 1)
  
  exp_WA <- exp(ww %*% alpha)
  pA_den <- opp %*% exp_WA
  # browser()
  # Part of MH ratio from P(A|O,alpha)
  pA_denstar <- pA_den +
    colSums(matrix(exp_WA[new] * plusminus, ncol = n_i)) # recall: matrix() lays out by column
  # Part of MH ratio from P(O|beta)  (logistic model)
  xb <- (xx %*% beta)[indnew]
  
  return(log(pA_den) - log(pA_denstar) + colSums(matrix(plusminus * xb, ncol = n_i)))
}

#' log_mh alpha
logmh_alpha <- function(alpha, alphastar, ww, opp, wa, prior) {
  
  # Calculate likelihood ratio
  exp_WA <- exp(ww %*% alpha)
  pA_den <- opp %*% exp_WA
  exp_WA_star <-  exp(ww %*% alphastar)
  pA_denstar <- opp %*% exp_WA_star
  
  logmh_alpha <- sum(wa * (alphastar - alpha)) + sum(log(pA_den) - log(pA_denstar)) +
    mvtnorm::dmvnorm(alphastar, prior$alpha$mu, solve(prior$alpha$Tau), log = TRUE) -
    mvtnorm::dmvnorm(alpha, prior$alpha$mu, solve(prior$alpha$Tau), log = TRUE)
  return(logmh_alpha)
}

#' log_mh beta
logmh_beta <- function(beta, betastar, xx, opp, mu_beta, Tau_beta) {
  XB <- xx %*% beta
  XB_star <- xx %*% betastar
  # Calculate likelihood ratio (using logistic structure)(don't count unemp.
  lrat <- sum((opp * (XB_star - XB)))    # the `canonical' part (unemp. cancels)
  logmh_beta <- lrat + sum(log(1 + exp(XB)) - log(1 + exp(XB_star))) +
    sum(mvtnorm::dmvnorm(t(betastar), mu_beta, solve(Tau_beta), log = TRUE)) -
    sum(mvtnorm::dmvnorm(t(beta), mu_beta, solve(Tau_beta), log = TRUE))
  return(logmh_beta)
}

cond_mu_beta <- function(Tau_beta, prior, beta) {
  n_j <- ncol(beta)
  Vinv <- prior$mu_beta$Tau + n_j * Tau_beta
  V <- solve(Vinv)
  m <- c(V %*% (prior$mu_beta$Tau %*% prior$mu_beta$mu +
                  n_j * Tau_beta %*% apply(beta, 1, mean)))
  return(list(m=m, V=V))
}

cond_Tau_beta <- function(mu_beta, prior, beta) {
  n_j <- ncol(beta)
  nu <- prior$Tau_beta$nu + n_j
  S <- solve(prior$Tau_beta$Sinv) + (beta - mu_beta) %*% t(beta - mu_beta)
  return(list(nu = nu, Sinv = solve(S)))
}

# ---- Adaptive MCMC utilities ----

#' Calculate the vcov matrix recursively
#' @param X_sample the latest sample of X
calculate_C <- function(X_sample, t, C, Xbar, sd, eps = 0.01) {
  d <- length(X_sample)
  X_sample <- c(X_sample) # Turn 1-row matrix into a vector
  
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

# ---- MCMC runs ----

match2sided <- function(iter, t0 = iter / 10, thin = 10,
                        C_alpha, C_beta,
                        starting_alpha, starting_beta,
                        frac_opp, prior,
                        ww, xx, choice, opp,
                        to_save = c("alpha", "beta"),
                        file, write = FALSE) {
  if (write & (t0 < iter)) {
    stop("Can't write intermediate result and do adaptive MCMC for now")
  }
  
  n_i <- nrow(xx)
  n_j <- nrow(ww)
  p_i <- ncol(xx)
  p_j <- ncol(ww)
  eps <- 0.01 # Small constant for beta proposal
  
  if (missing(prior)) {
    prior <- list(alpha = list(mu = rnorm(p_j),
                               Tau = solve(diag(abs(rnorm(p_j))))),
                  mu_beta = list(mu = rnorm(p_i),
                                 Tau = solve(diag(abs(rnorm(p_i))))),
                  Tau_beta = list(nu = p_i + 2,
                                  Sinv = solve(diag(abs(rnorm(p_i))))))
  }
  
  # ---- Starting values ----
  if (missing(starting_alpha)) {
    alpha <- rep(0, p_j)            # worker preferences  
  } else {
    alpha <- starting_alpha
  }
  alphastar <- alpha
  exp_WA <- exp(ww %*% alpha) # linear predictors of
  pA_den <- opp %*% exp_WA # vector of denominators in p(A | O, alpha)
  
  # beta starting values (from 1-sided logit estimates)
  if (missing(starting_beta)) {
    beta <- matrix(0, p_i, n_j)   # employer preferences
    for(j in 2:n_j) {
      y <- as.numeric(opp[, j])
      mod <- glm(y ~ . - 1, family=binomial,
                 data=as.data.frame(xx) )
      beta[, j] <- mod$coef
    }  
  } else {
    beta <- starting_beta
  }
  betastar <- beta
  
  XB <- xx %*% beta # worker side linear predictors (big matrix), n_i x n_j, same as opp
  
  # mu, Sigma starting values
  mu_beta <- prior$mu_beta$mu
  Tau_beta <- prior$Tau_beta$Sinv
  
  # ---- Pre compute ----
  tmp <- as.matrix(ww[choice, ])
  wa <- apply(tmp, 2, sum) # sum of characteristics of accepted jobs; used in alpha update
  
  sd_alpha <- (2.4 ** 2) / p_j / 10 # Scaling factor for alpha proposal
  sd_beta <- (2.4 ** 2) / (p_i * (n_j - 1) ) / 1000 # Scaling factor for beta proposal, minus 1 for unemployment
  
  # ---- Initialize storage ----
  acceptance_rate <- rep(0, 3)                  # Metropolis acceptance rates
  if ("opp" %in% to_save) {
    oppsave <- array(NA, dim = c(iter, n_i, n_j))
  } else {
    oppsave <- NA
  }
  asave <- matrix(NA, iter, p_j)   # saved alphas
  astarsave <- matrix(NA, iter, p_j)
  bsave <- array(NA, dim = c(iter, p_i, n_j)) # saved betas
  bstarsave <- array(NA, dim = c(iter, p_i, n_j)) # saved betas proposal
  mu_betasave <- matrix(NA, iter, p_i)
  Tau_betasave <- array(NA, dim = c(iter, p_i, p_i))
  logpost <- matrix(NA, iter, 3) # posterior density, i.e. joint, lp_A, lp_O
  ok <- matrix(NA, iter, 3) # ok_opp, ok_alpha, ok_beta
  
  # Naming the storage
  colnames(logpost) <- c("lpost", "lp_A", "lp_O")
  colnames(ok) <- c("mean(opp)", "alpha", "beta")
  if (!is.null(colnames(ww))) colnames(asave) <- colnames(ww)
  if (!is.null(colnames(ww))) colnames(astarsave) <- colnames(ww)
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
  
  # Store results
  logpost[1, 1] <- joint_lpdf(opp, choice,
                              alpha, beta, mu_beta, Tau_beta, prior,
                              ww, xx)
  logpost[1, 2] <- f_logp_A(opp, choice, alpha, ww)
  logpost[1, 3] <- f_logp_O(opp, beta, xx)
  ok[1, ] <- c(0, 0, 0)
  if ("opp" %in% to_save) oppsave[1, , ] <- opp
  asave[1, ] <- alpha
  astarsave[1, ] <- alphastar
  bsave[1, , ] <- beta # vectorize
  bstarsave[1, , ] <- betastar
  mu_betasave[1, ] <- mu_beta
  Tau_betasave[1, ,] <- Tau_beta

  # ---- Loop ----
  start <- Sys.time()
  for (i in 2:iter) {
    for (i_thin in 1:thin) {
      # Update opp
      num_new_offers <- floor(frac_opp * n_j) # per i
      new <- new_offer(n_i, n_j, num_new_offers)
      ind <- cbind(rep(1:n_i, each = num_new_offers), new)
      
      my_logmh_opp <- logmh_opp(opp, new, alpha, beta, ww, xx)
      ok_opp <- log(runif(n_i)) <= my_logmh_opp
      
      # Don't change an offer for an accepted job
      accepted_jobs_sampled <- colSums(matrix(new == rep(choice, each = num_new_offers),
                                              nrow = num_new_offers)) > 0
      ok_opp[accepted_jobs_sampled] <- F  # don't change an offer for an accepted job
      if (any(ok_opp)) {
        opp[ind][rep(ok_opp, each = num_new_offers)] <- !(opp[ind][rep(ok_opp, each=num_new_offers)]) # Update the opportunity set
        acceptance_rate[1] <- acceptance_rate[1] + mean(ok_opp)
      }
      
      # Update alpha
      if (i <= t0) {
        C_alpha_est <- C_alpha
      } else if (i > t0) {
        if (i == t0 + 1) {
          browser()
          # Calculate Cs and Xbars for the first time
          idx_first_accept_alpha <- min(which(ok[, "alpha"] == 1))
          alpha_samples <- asave[idx_first_accept_alpha:(i - 1), , drop = FALSE]
          C_alpha_est <- sd_alpha * var(alpha_samples) + sd_alpha * eps * diag(p_j)
          Xbar_alpha_est <- colMeans(alpha_samples)
        } else {
          if (i %% 200 == 0) { # Periodically check acceptance rate
            ok_rate_alpha <- mean(ok[(i - 200):(i - 1), "alpha"])
            if (ok_rate_alpha < 0.1) {
              sd_alpha <- sd_alpha * 1 / 2
              cat("sd_alpha decreases to", sd_alpha, "\n")
            } else if (ok_rate_alpha > 0.4) {
              sd_alpha <- sd_alpha * 2
              cat("sd_alpha increases to", sd_alpha, "\n")
            }
          }
          # Calculate Cs and Xbars recursively
          # (notice we're passing in old values of Cs and Xbars)
          res <- calculate_C(X_sample = asave[i - 1, , drop = FALSE], t = i, 
                             C = C_alpha_est, Xbar = Xbar_alpha_est,
                             sd = sd_alpha)
          Xbar_alpha_est <- res$Xbar
          C_alpha_est <- res$C
        }
      }
      alphastar <- alpha + c(rmvnorm(1, sigma = C_alpha_est))
      
      my_logmh_alpha <- logmh_alphaC(alpha, alphastar, ww, opp, wa,
                                     prior$alpha$mu, prior$alpha$Tau)
      ok_alpha <- ifelse(log(runif(1)) <= my_logmh_alpha, T, F)
      if (ok_alpha) {
        alpha <- alphastar
        acceptance_rate[2] <- acceptance_rate[2] + 1
      }
      
      # Update beta (except the beta of unemployment)
      if (i <= t0) {
        C_beta_est <- as.matrix(Matrix::bdiag(rep(list(C_beta), n_j - 1)))
      } else if (i > t0) {
        if (i == t0 + 1) {
          # Calculate Cs and Xbars for the first time
          idx_first_accept_beta <- min(which(ok[, "beta"] == 1))
          beta_samples <- bsave[idx_first_accept_beta:(i - 1), , 2:n_j, drop = FALSE]
          dim(beta_samples) <- c(i - idx_first_accept_beta, p_i * (n_j - 1))
          C_beta_est <- sd_beta * var(beta_samples) + sd_beta * eps * diag(p_i * (n_j - 1))
          Xbar_beta_est <- colMeans(beta_samples)
        } else {
          if (i %% 200 == 0) { # Periodically check acceptance rate
            ok_rate_beta <- mean(ok[(i - 200):(i - 1), "beta"])
            if (ok_rate_beta < 0.1) {
              sd_beta <- sd_beta * 1 / 10
              cat("sd_beta decreases to", sd_beta, "\n")
            } else if (ok_rate_beta > 0.4) {
              sd_beta <- sd_beta * 10
              cat("sd_beta increases to", sd_beta, "\n")
            }
          }
          # Calculate Cs and Xbars recursively
          # (notice we're passing in old values of Cs and Xbars)
          beta_sample <- c(bsave[i - 1, , 2:n_j, drop = FALSE])
          res <- calculate_C(X_sample = beta_sample, t = i, 
                             C = C_beta_est, Xbar = Xbar_beta_est,
                             sd = sd_beta)
          Xbar_beta_est <- res$Xbar
          C_beta_est <- res$C
        }
      }
      deviation <- matrix(rmvnorm(1, sigma = C_beta_est), nrow = p_i, ncol = (n_j - 1))
      deviation <- cbind(rep(0, p_i), deviation) # add the deviation for unemployment, which is 0
      betastar <- beta + deviation
      my_logmh_beta <- logmh_betaC(beta, betastar, xx, opp,
                                   mu_beta, Tau_beta)
      ok_beta <- ifelse(log(runif(1)) <= my_logmh_beta, T, F)
      if (ok_beta) {
        beta <- betastar
        acceptance_rate[3] <- acceptance_rate[3] + 1
      }
      
      # Update mu (multivariate normal)
      mu_beta_posterior <- cond_mu_beta(Tau_beta, prior, beta)
      mu_beta <- mvtnorm::rmvnorm(1, mu_beta_posterior$m,
                                  mu_beta_posterior$V)
      mu_beta <- c(mu_beta) # Flatten into a vector
      
      # Update Tau (Wishart)
      Tau_beta_posterior <- cond_Tau_beta(mu_beta, prior, beta)
      Tau_beta <- MCMCpack::rwish(Tau_beta_posterior$nu,
                                  Tau_beta_posterior$Sinv)
      
    }
    
    # Store results
    
    logpost[i, 1] <- joint_lpdf(opp, choice,
                                alpha, beta, mu_beta, Tau_beta, prior,
                                ww, xx)
    logpost[i, 2] <- f_logp_A(opp, choice, alpha, ww)
    logpost[i, 3] <- f_logp_O(opp, beta, xx)
    ok[i, ] <- c(mean(ok_opp), ok_alpha, ok_beta)
    if ("opp" %in% to_save) oppsave[i, , ] <- opp
    asave[i, ] <- alpha
    astarsave[i, ] <- alphastar
    bsave[i, , ] <- beta # vectorize
    bstarsave[i, , ] <- betastar
    mu_betasave[i, ] <- mu_beta
    Tau_betasave[i, ,] <- Tau_beta  
    
    # Increment
    interval_length <- 1000
    if (i %% interval_length == 0) {
      cat("Iteration", i, "done", Sys.time() - start, "\n")
      start <- Sys.time()
      if (write) {
        write_to_disk(idx = i, interval_length = interval_length,
                      vars_to_write = c("asave" = "alpha", "astarsave" = "alphastar",
                                        "bsave" = "beta", "bstarsave" = "betastar"))
      }
    }
  }
  
  return(list(opp = oppsave, alpha = asave, beta = bsave,
              alphastar = astarsave, betastar = bstarsave,
              mu_beta = mu_betasave, Tau_beta = Tau_betasave,
              lp = logpost, ok = ok,
              acceptance_rate = acceptance_rate / iter,
              mcmc_settings = list(iter = iter, frac_opp = frac_opp, thin = thin,
                                   C_alpha = C_alpha, C_beta = C_beta)))
}