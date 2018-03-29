#' Density plot for posterior samples
#' @param data A matrix of n_sim x n_coef. colnames(data) will be used for label
plot_posterior_dens <- function(data, coefnames, xlab, order = FALSE) {
  
  dists <- vector("list", length = ncol(data))
  for (i in 1:ncol(data)) {
    dens <- density(data[, i])
    dists[[i]] <- data.frame(x = dens[['x']], y = dens[['y']] / max(dens[['y']]) * 0.75)
  }
  
  quantiles <- vector("list", length = ncol(data))
  for (i in 1:ncol(data)) {
    hpd <- coda::HPDinterval(coda::mcmc(data[, i]), prob = 0.95)
    dens <- density(data[, i], from = hpd[1], to = hpd[2])
    quantiles[[i]] <- data.frame(x = dens[['x']], y = dens[['y']] / max(dens[['y']]) * 0.75)
  }
  
  if (missing(coefnames)) {
    coefnames <- colnames(data)
  }
  
  if (order) {
    dists <- dists[order(colMeans(data))]
    quantiles <- quantiles[order(colMeans(data))]
    coefnames <- coefnames[order(colMeans(data))]
  }
  
  # blank canvas
  par(mar = c(5.1, 8.1, 0.1, 1.1))
  plot(0, 0, type="n",
       xlim=c(min(data), max(data)),
       ylim=c(1, ncol(data) + 1),
       xlab= ifelse(missing(xlab), "Coefficient", xlab),
       ylab="", yaxt="n")
  for (ii in 1:ncol(data)) {
    between_x <- quantiles[[ii]][ ,1]
    between_y <- c(0, quantiles[[ii]][ , 2], 0)
    between_x <- c(min(between_x), between_x, max(between_x))
    polygon(dists[[ii]][,1], dists[[ii]][,2]+ii, col="white", border=NA)
    polygon(between_x, between_y+ii, col="#CCCCCC", border=NA) # grey shade
    polygon(dists[[ii]][,1], dists[[ii]][,2]+ii) # density plot
  }
  abline(v=0, lty=2)
  axis(2, at = 1:ncol(data), coefnames, las=1) # y tick label
}

#' Simulate the observed choice, one simulation
#' @param alpha i's preferences
#' @param beta j's preferences
#' @return a n_i vector of the choices
f_sim <- function(alpha, beta) {
  # countries make offers to firms
  alpha1 <- alpha[sample(1:nrow(alpha), 1), ]
  beta1 <- beta[sample(1:nrow(beta), 1), , ]
  
  country_utilities <- (xx %*% beta1) + rlogis(n = n_i * n_j) # linear pred + error
  opp <- country_utilities > 0
  
  wa <- alpha1 %*% t(ww) # linear pred
  mnc_utilties <- t(replicate(n_i, wa + evd::rgumbel(n = n_j), simplify = TRUE))
  mnc_options <- opp * mnc_utilties
  mnc_choices <- rownames(ww)[apply(mnc_options, 1, which.max)]
  return(mnc_choices)
}

#' Visualize sim result
#' @param sim_result a matrix of n_sim x n_j, value being the test quantities
#' @param observed a data frame of the observed quantities for j
f_plot_sim <- function(sim_result, observed) {
  molten <- reshape2::melt(sim_result, variable.name = "choices") %>%
    group_by(choices) %>%
    summarise(mean = mean(value, na.rm = TRUE), 
              lower95 = quantile(value, 0.025, na.rm = TRUE),
              upper95 = quantile(value, 0.975, na.rm = TRUE))
  pd <- molten %>% inner_join(observed, by = "choices")
  ggplot(pd, aes(x = choices)) +
    geom_pointrange(aes(y = mean, ymin = lower95, ymax = upper95,
                        color = "Predicted")) +
    geom_point(aes(y = value, color = "Observed"))
}
  
#' Coeficient plot with pointrange, ordering the coef based on posterior mean
#' @param data n_sim x n_coef matrix
plot_posterior_coef <- function(data) {
  posterior_summary <- data.frame(mean = colMeans(data),
                                  coda::HPDinterval(mcmc(data), prob = 0.95))
  posterior_summary$coef <- row.names(posterior_summary)
  
  ggplot(data = posterior_summary,
         aes(fct_reorder(factor(coef), mean))) + # reorder according to mean posterior
    geom_pointrange(aes(y = mean, ymin = lower, ymax = upper)) +
    geom_hline(aes(yintercept = 0), linetype = 2) +
    coord_flip() + theme_bw()
}

#' Plot the predicted probability of being chosen, multinomial
#' @param mcmc n_sim x n_coef matrix
#' @param ww n_j x p_j matrix of covariates
#' @param observed a n_j x 2 data frame of observed proportion of the choices
#' must have column name "choices"
plot_multinom_pred <- function(alpha, ww, observed) {
  lin_pred <- alpha %*% t(ww) # linear predictors
  
  probs <- exp(lin_pred) / rowSums(exp(lin_pred))
  probs_molten <- reshape2::melt(as_tibble(probs),
                                 variable.name = "choices") %>%
    group_by(choices) %>%
    summarise(mean = mean(value), 
              lower95 = quantile(value, 0.025),
              upper95 = quantile(value, 0.975))
  pd <- probs_molten %>% inner_join(observed, by = "choices")
  ggplot(pd, aes(x = choices)) +
    geom_pointrange(aes(y = mean, ymin = lower95, ymax = upper95,
                        color = "Predicted")) +
    geom_point(aes(y = value, color = "Observed")) +
    labs(y = "Probability")
}
