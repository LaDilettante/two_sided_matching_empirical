rm(list = ls())
source("0_plot_functions.R")
library("tidyverse")

# ---- Load results ----

res <- readRDS("../result/labor_occ5_2018-04-03 13:17:02.RData")
start <- res$mcmc_settings$iter / 2
end <- res$mcmc_settings$iter
xx <- res$data[[1]]
ww <- res$data[[2]]
n_i <- nrow(res$data[[1]])
n_j <- nrow(res$data[[2]])

alpha <- res$alpha[start:end, ]
beta <- res$beta[start:end, , ]

# ---- Coefficient plot ----

pdf("../figure/labor_occ5_alpha.pdf", w = 4.5, h = 4.5)
plot_posterior_dens(alpha,
                    coefnames = c("Autonomy" = "aut", "Prestige" = "pres"))
dev.off()
colMeans(alpha)

pdf("../figure/labor_occ5_beta_educ_age.pdf", w = 7.5, h = 4.4)
par(mfrow = c(1, 2))
plot_posterior_dens(beta[, 'educ', 2:n_j], main = "Preference for education",
                    coefnames = c("Services" = "Sales, Clerical, Services"))
plot_posterior_dens(beta[, 'age', 2:n_j], main = "Preference for age",
                    coefnames = c("Services" = "Sales, Clerical, Services"))
par(mfrow = c(1, 1))
dev.off()

colMeans(res$beta[start:end, 'educ', ])
colMeans(res$beta[start:end, 'age', ])

# ---- Predicted probability of hiring ----

pd1 <- f_predicted_effect(varname = "educ", 
                            varvalues = quantile(xx[ , "educ"], probs = seq(0, 1, by = 0.05)),
                            actorj = "Professional")
pd2 <- f_predicted_effect(varname = "educ", 
                         varvalues = quantile(xx[ , "educ"], probs = seq(0, 1, by = 0.05)),
                         actorj = "Sales, Clerical, Services")
pd <- dplyr::bind_rows(pd1, pd2) %>%
  mutate(educ = educ + mean(df_employee_raw$educ)) # add back the mean

ggplot(pd, aes(x = educ)) +
  geom_line(aes(y = mean, color = actorj)) +
  geom_ribbon(aes(ymin = lower95, ymax = upper95, fill = actorj), 
              alpha = 0.25) +
  scale_color_discrete("") + scale_fill_discrete("") +
  labs(x = "Years of education", y = "Probability of hiring")
ggsave("../figure/labor_occ5_educ_effect_on_hiring.pdf", w = 7, h = 3.5)

# ---- Convergence ----

pdf("../figure/labor_occ5_misc_traceplots.pdf", w = 7, h = 5.5)
par(mfrow = c(2, 2))
plot(res$lp[, 1], type='l', xlab = 'Iteration', ylab = '', main = 'posterior density')
plot(res$alpha[ , 'pres'], type = 'l', xlab = 'Iteration', ylab = '',
     main = "Workers' pref for prestige")
plot(res$beta[ , 'educ', 'MFG Blue Collar'], type = 'l', 
     xlab = 'Iteration', ylab = '', main = "MFG's pref for education")
plot(res$beta[ , 'nonwhite', 'MFG Blue Collar'], type = 'l', 
     xlab = 'Iteration', ylab = '', main = "MFG's pref for non-white")
par(mfrow = c(1, 1))
dev.off()

# ---- opp ----

tmp <- t(apply(res$opp, 1, function(mat) colMeans(mat)))
plot(tmp[, 1], type = 'l')
plot(tmp[, 2] - colMeans(obs_opp)[2], type = 'l')
plot(tmp[, 3] - colMeans(obs_opp)[3], type = 'l')
plot(tmp[, 4] - colMeans(obs_opp)[4], type = 'l')
