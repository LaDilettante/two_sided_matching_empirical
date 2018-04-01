# --- Require df_mnc ---

rm(list = ls())
source("0_plot_functions.R")
library("tidyverse")

# result <- readRDS("../result/japan96_2018-03-22 20_49_28.RData")
# result <- readRDS("../result/japan96_2018-03-21 21_36_08.RData")
# result <- readRDS("../result/japan96_2018-03-23 19_28_10.RData")
result <- readRDS("../result/japan96_2018-03-24 12_34_01.RData")
start <- result$mcmc_settings$iter / 2
end <- result$mcmc_settings$iter
xx <- result$data[[1]]
ww <- result$data[[2]]
n_i <- nrow(result$data[[1]])
n_j <- nrow(result$data[[2]])

alpha <- result$alpha[start:end, ]
beta <- result$beta[start:end, , ]

# ---- Coefficient plot ----

source("0_plot_functions.R")

pdf("../figure/japan96_alpha.pdf", w = 5, h = 5.5)
plot_posterior_dens(alpha, coefnames = c("log GDP", "log GDP per cap", "human capital\nindex"))
dev.off()

pdf("../figure/japan96_beta_ltemp_luscptl.pdf", w = 9, h = 5)
par(mfrow = c(1, 2))
plot_posterior_dens(result$beta[start:end, 'ltemp', ], 
                    xlab = "Coefficient for log Employment", order = TRUE)
plot_posterior_dens(result$beta[start:end, 'luscptl', ], 
                    xlab = "Coefficient for log Capital", order = TRUE)
dev.off()


pdf("../figure/japan96_beta_intrd_intexp.pdf", w = 9, h = 5)
par(mfrow = c(1, 2))
plot_posterior_dens(result$beta[start:end, 'int_r_d', ], 
                    xlab = "Coefficient for R&D intensity", order = TRUE)
plot_posterior_dens(result$beta[start:end, 'int_exp', ], 
                    xlab = "Coefficient for export intensity", order = TRUE)
dev.off()

# ---- Predicted effect of a variable on probability of a country offering to a firm ----

f_predicted_effect <- function(varname, varvalues, actorj) {
  beta_tmp <- beta[sample(1:nrow(beta), 1000, replace = TRUE), , actorj]
  
  scenario <- t(replicate(length(varvalues), apply(xx, 2, median)))
  scenario[ , varname] <- varvalues
  lin_pred <- beta_tmp %*% t(scenario)
  prob <- exp(lin_pred) / (1 + exp(lin_pred))
  molten <- reshape2::melt(prob)
  molten %>% group_by(Var2) %>% 
    summarise(mean = mean(value),
              lower95 = quantile(value, probs = 0.05),
              upper95 = quantile(value, probs = 0.95)) %>%
    mutate(actorj = actorj, !!varname := varvalues)
}

china <- f_predicted_effect(varname = "ltemp", 
                   varvalues = quantile(xx[ , "ltemp"], probs = seq(0, 1, by = 0.05)),
                   actorj = "China")
sk <- f_predicted_effect(varname = "ltemp", 
                   varvalues = quantile(xx[ , "ltemp"], probs = seq(0, 1, by = 0.05)),
                   actorj = "South Korea")
pd <- dplyr::bind_rows(china, sk) %>%
  mutate(temp = exp(ltemp))

pdf("../figure/japan96_effect_of_temp.pdf", w = 9, h = 5)
ggplot(pd, aes(x = ltemp)) +
  geom_line(aes(y = mean, color = actorj)) +
  geom_ribbon(aes(ymin = lower95, ymax = upper95, fill = actorj), alpha = 0.25) +
  geom_rug(data = df_mnc) +
  scale_x_continuous(breaks = log(c(12.5, 50, 200, 800, 3200)),
                     labels = function(x) exp(x)) +
  scale_color_discrete("Country") +
  scale_fill_discrete("Country") + 
  labs(x = "Number of Employees", y = "Probability")
dev.off()

ggplot(pd, aes(x = temp)) +
  geom_line(aes(y = mean, color = actorj)) +
  geom_ribbon(aes(ymin = lower95, ymax = upper95, fill = actorj), alpha = 0.25) +
  geom_rug(data = df_mnc, aes(x = exp(ltemp)))

p1 <- ggplot(mpg, aes(displ, hwy)) +
  geom_point()
p1
p1 + scale_y_log10()

# ---- prob of a country being chosen, predicted vs true ----

# but this doesn't take into account the other side though
observed_prop <- df_mnc %>% dplyr::count(nation) %>% 
  mutate(value = n / sum(n)) %>%
  rename(choices = nation)

pdf("../figure/japan96_prob_being_chosen_by_MNC_one_sided.pdf", w = 5.5, h = 3)
plot_multinom_pred(alpha, ww, observed_prop) +
  labs(x = "Country") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()

# ---- prob a country being chosen, simulate both sides ----
sim_choice <- t(replicate(1000, f_sim(alpha, beta)))
# Coerce to data frame so that we can use dplyr::bind_rows,
# which fills missing with NA
probs_country <- apply(sim_choice, 1, 
                       function(x) as.data.frame(as.list(table(x) / n_i))) %>%
  dplyr::bind_rows()
pdf("../figure/japan96_prob_being_chosen_by_MNC_two_sided.pdf", w = 5.5, h = 3)
f_plot_sim(probs_country, observed_prop) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(x = "Country", y = "Probability")
dev.off()

# ---- Profile of firms in a country ----

f_profile <- function(choice) {
  pd <- cbind.data.frame(df_mnc, sim_choice = choice)
  tmp <- pd %>% group_by(sim_choice) %>%
    summarise(mean(ltemp), mean(luscptl), mean(int_r_d), mean(int_exp),
              median(ltemp), median(luscptl), median(int_r_d), median(int_exp),
              var(ltemp), var(luscptl), var(int_r_d), var(int_exp))
  ret <- tmp %>% select(-sim_choice) %>% as.matrix()
  rownames(ret) <- tmp$sim_choice
  return(ret)
}
sim_profile <- apply(sim_choice, 1, f_profile)

# Mean of ltemp
sim_result <- dplyr::rbind_list(lapply(sim_profile, function(mat) mat[ , "mean(ltemp)"]))
observed_mean_ltemp <- df_mnc %>% group_by(nation) %>% 
  summarize(value = mean(ltemp)) %>%
  rename(choices = nation)

pdf("../figure/japan96_sim_mean_ltemp.pdf", w = 5.5, h = 3)
f_plot_sim(sim_result, observed_mean_ltemp) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(x = "Country", y = "Mean of log Employee")
dev.off()

# Variance of ltemp
sim_result <- dplyr::rbind_list(lapply(sim_profile, function(mat) mat[ , "var(ltemp)"]))
observed_var_ltemp <- df_mnc %>% group_by(nation) %>% 
  summarize(value = var(ltemp)) %>%
  rename(choices = nation)
pdf("../figure/japan96_sim_var_ltemp.pdf", w = 5.5, h = 3)
f_plot_sim(sim_result, observed_var_ltemp) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(x = "Country", y = "Variance of log Employee")
dev.off()


# ---- Predicted location of a particular firm ----

f_sim_firm_location <- function() {
  alpha1 <- alpha[sample(1:nrow(alpha), 1), ]
  beta1 <- beta[sample(1:nrow(beta), 1), , ]
  
  country_utilities <- (xx %*% beta1) + rlogis(n = n_i * n_j) # linear pred + error
  opp <- country_utilities > 0
  
  wa <- alpha1 %*% t(ww) # linear pred
  mnc_utilties <- t(replicate(n_i, wa + evd::rgumbel(n = n_j), simplify = TRUE))
  mnc_options <- opp * mnc_utilties
  mnc_choices <- rownames(ww)[apply(mnc_options, 1, which.max)]
  df_mnc$nation == mnc_choices
}

sim_result <- t(replicate(1000, f_sim_firm_location()))


# ---- Predicted effect of hc (hold median for other variables) ----

# probability of choosing a country based on different levels of hc

# let's take Thailand, vary its level of hc and 
# show how the probability of choosing Thailand differs

set up scenario
calculate probs of being chosens

ww_scenario <- ww
ww_scenario['Thailand', 'hc'] <- 2.6
sims <- alpha[sample(1:nrow(alpha), replace = TRUE, size = 1000), ]
lin_pred <- sims %*% t(ww_scenario) # linear predictors
probs <- exp(lin_pred) / rowSums(exp(lin_pred))


# ---- Posterior predictive check: similarity of firms in a country ----



# heatmap of opp changing

# 
