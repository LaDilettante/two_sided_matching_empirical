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
  # browser()
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

# ---- MCMC runs ----

match2sided <- function(iter, eps_alpha, eps_beta, frac_beta, frac_opp,
                        ww, xx, choice, opp) {
  n_i <- nrow(xx)
  n_j <- nrow(ww)
  p_i <- ncol(xx)
  p_j <- ncol(ww)
  
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
  beta  <- matrix(0, p_i, n_j)   # employer preferences
  for(j in 1:n_j) {
    y <- as.numeric(opp[, j])
    mod <- glm(y ~ . - 1, family=binomial,
               data=as.data.frame(xx) )
    beta[, j] <- mod$coef
  }
  XB <- xx %*% beta # worker side linear predictors (big matrix), n_i x n_j, same as opp
  
  # mu, Sigma starting values
  mu_beta <- prior$mu_beta$mu
  Tau_beta <- prior$Tau_beta$Sinv
  
  # ---- Pre compute ----
  tmp <- as.matrix(ww[choice, ])
  wa <- apply(tmp, 2, sum) # sum of characteristics of accepted jobs; used in alpha update
  
  bmat <- eps_beta * matrix(1, p_i, n_j)
  
  # ---- Initialize storage ----
  acceptance_rate <- c("opp" = 0, "alpha" = 0, "beta" = 0) # Metropolis acceptance rates
  asave <- matrix(NA, iter, p_j)   # saved alphas
  bsave <- array(NA, dim = c(iter, p_i, n_j)) # saved betas
  mu_betasave <- matrix(NA, iter, p_i)
  Tau_betasave <- array(NA, dim = c(iter, p_i, p_i))
  logpost <- matrix(NA, iter, 3) # posterior density, i.e. lp_A, lp_O, joint
  
  # ---- Loop ----
  for (i in 1:iter) {
    # Update opp
    # browser()
    num_new_offers <- floor(frac_opp * n_j) # per i
    new <- replicate(n_i,
                     sample(2:n_j, size=floor(frac_opp * n_j), replace = FALSE))
    new <- c(new) # Flatten 1-column matrix into a vector
    ind <- cbind(rep(1:n_i, each = num_new_offers), new)
    
    my_logmh_opp <- logmh_opp(opp, new, alpha, beta, ww, xx)
    ok_opp <- log(runif(n_i)) <= my_logmh_opp
    
    # Don't change an offer for an accepted job
    accepted_jobs_sampled <- colSums(matrix(new == rep(choice, each = num_new_offers),
                                            nrow = num_new_offers)) > 0
    ok_opp[accepted_jobs_sampled] <- F  # don't change an offer for an accepted job
    if (any(ok_opp)) {
      opp[ind][rep(ok_opp, each = num_new_offers)] <- !(opp[ind][rep(ok_opp, each=num_new_offers)]) # Update the opportunity set
      acceptance_rate[1] <- acceptance_rate[1] + 1
    }
    
    # Update alpha
    deviation <- eps_alpha * runif(p_j, min=-1, max=1) # Symmetric proposal
    alphastar <- alpha + deviation
    
    my_logmh_alpha <- logmh_alpha(alpha, alphastar, ww, opp, wa, prior)
    ok_alpha <- ifelse(log(runif(1)) <= my_logmh_alpha, T, F)
    if (ok_alpha) {
      alpha <- alphastar
      acceptance_rate[2] <- acceptance_rate[2] + 1
    }
    
    # Update beta
    whichones <- sample(c(0, 1), size = p_i * n_j, replace = TRUE,
                        prob = c(1 - frac_beta, frac_beta))
    # Sample betastar from a [-eps_beta, eps_beta] box around beta
    rmat <- matrix(runif(p_i * n_j, min=-1, max=1) * whichones,
                   nrow = p_i, ncol = n_j)
    deviation <- bmat * rmat
    betastar <- beta + deviation
    my_logmh_beta <- logmh_beta(beta, betastar, xx, opp, mu_beta, Tau_beta)
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
    
    # Store results
    logpost[i, 1] <- joint_lpdf(opp, choice,
                              alpha, beta, mu_beta, Tau_beta, prior,
                              ww, xx)
    logpost[i, 2] <- f_logp_A(opp, choice, alpha, ww)
    logpost[i, 3] <- f_logp_O(opp, beta, xx)
    asave[i, ] <- alpha
    bsave[i, ,] <- beta # vectorize
    mu_betasave[i, ] <- mu_beta
    Tau_betasave[i, ,] <- Tau_beta
    
    # Increment
    if (i %% 100 == 0) cat("Iteration", i, "done", "\n")
  }
  
  if (!is.null(colnames(ww))) colnames(asave) <- colnames(ww)
  if (!is.null(colnames(xx))) {
    dimnames(bsave)[[2]] <- colnames(xx)
    colnames(mu_betasave) <- colnames(xx)
    dimnames(Tau_betasave)[[2]] <- colnames(xx)
    dimnames(Tau_betasave)[[3]] <- rownames(xx)
  }
  if (!is.null(rownames(ww))) {
    dimnames(bsave)[[3]] <- rownames(ww)
  }
  
  return(list(alpha = asave, beta = bsave,
              mu_beta = mu_betasave, Tau_beta = Tau_betasave,
              lp = logpost,
              acceptance_rate = acceptance_rate / iter,
              mcmc_settings = list(iter = iter,
                                   eps_alpha = eps_alpha,
                                   eps_beta = eps_beta,
                                   frac_beta = frac_beta,
                                   frac_opp = frac_opp)))
}