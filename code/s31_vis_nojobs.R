rm(list = ls())
source("0_plot_functions.R")
library("tidyverse")
library("coda")
library("stargazer")

res <- readRDS("../result/simlabor_nojobs_2018-03-31 19_02_45.RData")
start <- res$mcmc_settings$iter / 2
end <- res$mcmc_settings$iter
xx <- res$data[[1]]
ww <- res$data[[2]]
n_i <- nrow(res$data[[1]])
n_j <- nrow(res$data[[2]])

alpha <- res$alpha[start:end, ]
beta <- res$beta[start:end, , ]

# ----- Summary statistics -----

stargazer(py$df_employee %>% select(educ, age, nonwhite), 
          summary = TRUE,
          covariate.labels = c("Years of education", "Age", "Non-white"),
          float = FALSE,
          out = "../table/labor_occ5_summary_employee.tex")

stargazer(py$df_employer %>% select(occ5, pres, supervisor, aut) %>%
            mutate(occ5 = recode(occ5, 
                                 "MFG Blue Collar" = "Manufacturing Blue Collar")), 
          covariate.labels = c("Firm category", "Prestige", "Pr(Supervisor)", "Autonomy"),
          float = FALSE, summary = FALSE, rownames = FALSE,
          out = "../table/labor_occ5_summary_employer.tex")

# ---- Trace plot  ----

source("0_plot_functions.R")

# Trace plot alpha
pdf("../figure/sim_labor_nojobs_nounemp_alpha.pdf", w = 7, h = 4.5)
par(mfrow = c(1, 2))
traceplot(mcmc(res$alpha)[, 1], xlab = "Workers' preference for job's prestige")
abline(h = py$true_alpha[1], col = "red")
traceplot(mcmc(res$alpha)[, 2], xlab = "Workers' preference for job's autonomy")
abline(h = py$true_alpha[2], col = "red")
dev.off()

# Trace plot beta
for (j in 1:py$n_j) {
  pdf(paste0("../figure/sim_labor_nojobs_beta_emp", j, ".pdf"), w = 6.5, h = 5)
  par(mfrow = c(2, 2), mar = c(4, 4, 2.5, 1))
  xlabs <- c("Intercept", "Preference for education", 
             "Preference for age", "Preference for non-white")
  for (i in 1:4) {
    traceplot(mcmc(res$beta[, i, j]), xlab = xlabs[i])
    abline(h = py$true_beta_std[j, i], col = "red")
  }
  dev.off()
}

# ---- TSL vs CL ----

# Conditional logit
log_lik <- function(alpha) {
  wa <- py$ww %*% alpha
  lse_wa <- logsumexpC(wa)
  sum(wa[py$choice] - lse_wa)
}
m_multinom <- optim(c(0, 0), log_lik, 
                    control = list("fnscale" = -1), hessian = TRUE)
alpha_cl <- mvtnorm::rmvnorm(nrow(alpha), mean = m_multinom$par,
                              sigma = solve(-m_multinom$hessian))

# Poisson (identical result as CL)
df_count <- data.frame(table(py$choice),
                         py$df_ww)
m_pois <- summary(glm(Freq ~ pres + aut, family = "poisson", data = df_count))

# Negative binomial
m_nb <- MASS::glm.nb(Freq ~ pres + aut, data = df_count)

pd <- data.frame(alpha, alpha_cl) %>%
  rename(prestige_tsl = X1, autonomy_tsl = X2, 
         prestige_cl = X1.1, autonomy_cl = X2.1) %>% 
  gather() %>% separate(key, sep = "_", into = c("variable", "model"))

pdf("../figure/sim_labor_nojobs_alpha_tsl_vs_cl.pdf", w = 8, h = 3.5)
ggplot(pd) + geom_density(aes(x = value, fill = model), alpha = 0.5) +
  geom_vline(data = data.frame(variable = c("prestige", "autonomy"),
                               true_value = py$true_alpha),
             aes(xintercept = true_value, color = "true value")) +
  facet_wrap(~ variable, scales = "free") +
  scale_color_grey("") +
  scale_fill_hue("", labels = c("tsl" = "Two-sided logit", 
                                "cl" = "One-sided\nconditional logit")) +
  labs(x = "Coefficient")
dev.off()

# True opp vs obs choice
p_trueopp <- f_plot_opp(py$true_opp) +
  scale_x_continuous(breaks = 1:ncol(py$true_opp),
                     labels = py$df_employer$occ5) +
  theme(axis.text.x = element_text(angle = 15, vjust = 1, hjust=0.5)) +
  labs(x = "Firm index", y = "Worker index", title = "True opportunity set")

choice_matrix <- matrix(FALSE, nrow = py$n_i, ncol = py$n_j)
choice_matrix[cbind(1:py$n_i, py$choice)] <- TRUE
p_obschoice <- f_plot_opp(choice_matrix) +
  scale_x_continuous(breaks = 1:ncol(py$true_opp),
                     labels = py$df_employer$occ5) +
  theme(axis.text.x = element_text(angle = 15, vjust = 1, hjust=0.5)) +
  labs(x = "Firm index", y = "Worker index", title = "Observed choice")

ggpubr::ggarrange(p_trueopp, p_obschoice, ncol = 2, common.legend = TRUE) 
ggsave("../figure/sim_labor_nojobs_trueopp_obschoice.pdf", 
       width = 10.4, height = 6.8)

# ---- Visualize the opp ----

offer_rate <- t(apply(res$opp, 1, colMeans))
intercept <- res$beta[ , 1, ] # intercept

pdf("../figure/sim_labor_nojobs_nounemp_managerial.pdf", w = 8, h = 3.2)
par(mfrow = c(1, 3))
j <- 2 # Managerial
plot(offer_rate[, j], type = 'l', xlab = "Iterations", ylab = "Offer rate")
abline(h = colMeans(py$true_opp)[j], col = 'red')
traceplot(mcmc(res$beta[, 1, j]), ylab = "Intercept", ylim = c(-1.5, 4))
abline(h = py$true_beta_std[j, 1], col = "red")
traceplot(mcmc(res$beta[, 3, j]), ylab = "Preference for educ")
abline(h = py$true_beta_std[j, 3], col = "red")
dev.off()

my.plot.mcmc(mcmc(res$beta[, , 2]), parameters = py$true_beta_std[2, ])


pdf("../figure/sim_labor_nojobs_nounemp_offer_rate.pdf", w = 6.8, h = 3.8)
par(mfrow = c(1, 2))
plot(offer_rate[, 1], type = 'l', 
     xlab = "Iterations", ylab = "Offer rate", main = "Professional")
abline(h = colMeans(py$true_opp)[1], col = 'red')
plot(offer_rate[, 4], type = 'l', 
     xlab = "Iterations", ylab = "Offer rate", main = "Other blue collar")
abline(h = colMeans(py$true_opp)[4], col = 'red')
dev.off()

# ---- Correlation between beta and the opp ----

offer_rate <- t(apply(res$opp, 1, colMeans))
intercept <- res$beta[ , 1, ] # intercept

idx_highed <- which(py$xx[, 2] > quantile(py$xx[, 2], probs = c(0.75)))
offer_rate_higheduc <- t(apply(res$opp[, idx_highed, ], 1, colMeans))
beta_educ <- res$beta[, 2, ]

# Correlation of opp and beta for professional
start <- res$mcmc_settings$iter / 2
end <- res$mcmc_settings$iter
pdf("../figure/sim_labor_nojobs_opp_beta_correlation_managerial.pdf", 
    w = 8, h = 4.5)
par(mfrow = c(1, 2))
j <- 3 # Managerial
plot(offer_rate[start:end, j], intercept[start:end, j],
     cex = 0.5, col=rgb(0, 0, 0, 0.25),
     xlab = "Offer rate", ylab = "Intercept")
plot(offer_rate_higheduc[start:end, j], beta_educ[start:end, j],
     cex = 0.5, col=rgb(0, 0, 0, 0.25),
     xlab = "Offer rate for educated workers", ylab = "Preference for education")
dev.off()
