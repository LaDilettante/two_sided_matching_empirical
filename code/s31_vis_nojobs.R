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
pdf("../figure/sim_labor_nojobs_alpha.pdf", w = 7, h = 4.5)
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

log_lik <- function(alpha) {
  wa <- py$ww %*% alpha
  lse_wa <- logsumexpC(wa)
  sum(wa[py$choice] - lse_wa)
}
m_multinom <- optim(c(0, 0), log_lik, 
                    control = list("fnscale" = -1), hessian = TRUE)
alpha_cl <- mvtnorm::rmvnorm(nrow(alpha), mean = m_multinom$par,
                              sigma = solve(-m_multinom$hessian))

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