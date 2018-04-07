# --- Require df_mnc ---

# rm(list = ls())
source("0_plot_functions.R")
library("tidyverse")

# res <- readRDS("../result/japan96_2018-04-06 15:24:11.RData")
res <- readRDS("../result/japan96_2018-04-05 00:25:23.RData")
start <- res$mcmc_settings$iter / 2
end <- res$mcmc_settings$iter
xx <- res$data[[1]]
xx_mean <- colMeans(df_mnc %>% select(ltemp, luscptl))
ww <- res$data[[2]]
n_i <- nrow(res$data[[1]])
n_j <- nrow(res$data[[2]])

alpha <- res$alpha[start:end, ]
alphaest <- colMeans(alpha)
beta <- res$beta[start:end, , ]

alphaest %*% t(res$data$ww) # See which country is most desirable

# ---- Coefficient plot ----

lookup_alpha_coef <- c("log GDP" = "lgdp", "log GDP per cap" = "lgdppc", 
                       "Human capital" = "hc")
lookup_beta_coef <- c("log Employment" = "ltemp", "log Capital" = "luscptl",
                      "R&D intensity" = "int_r_d", "Export intensity" = "int_exp")
  
pdf("../figure/japan96_alpha.pdf", w = 5, h = 5.5)
plot_posterior_dens(alpha, 
                    coefnames = lookup_alpha_coef)
dev.off()

pdf("../figure/japan96_beta_Taiwan_Indonesia.pdf", w = 7.4, h = 4.5)
par(mfrow = c(1, 2))
plot_posterior_dens(res$beta[start:end, 2:5, c("Taiwan")],
                    coefnames = lookup_beta_coef,
                    main = "Taiwan preference")
plot_posterior_dens(res$beta[start:end, 2:5, c("Indonesia")],
                    coefnames = lookup_beta_coef,
                    main = "Indonesia preference")
dev.off()

pdf("../figure/japan96_beta_ltemp_luscptl.pdf", w = 9, h = 5)
par(mfrow = c(1, 2))
plot_posterior_dens(res$beta[start:end, 'ltemp', ], 
                    xlab = "Coefficient for log Employment", order = TRUE)
plot_posterior_dens(res$beta[start:end, 'luscptl', ], 
                    xlab = "Coefficient for log Capital", order = TRUE)
dev.off()

pdf("../figure/japan96_beta_intrd_intexp.pdf", w = 9, h = 5)
par(mfrow = c(1, 2))
plot_posterior_dens(res$beta[start:end, 'int_r_d', ], 
                    xlab = "Coefficient for R&D intensity", order = TRUE)
plot_posterior_dens(res$beta[start:end, 'int_exp', ], 
                    xlab = "Coefficient for export intensity", order = TRUE)
dev.off()

# ---- Predicted effect of a variable on probability of a country offering to a firm ----

country1 <- f_predicted_effect(varname = "int_exp", 
                   varvalues = quantile(xx[ , "int_exp"], probs = seq(0, 1, by = 0.05)),
                   actorj = "Malaysia")
country2 <- f_predicted_effect(varname = "int_exp", 
                   varvalues = quantile(xx[ , "int_exp"], probs = seq(0, 1, by = 0.05)),
                   actorj = "Taiwan")
pd <- dplyr::bind_rows(country1, country2) %>%
  # add back the mean and scale
  mutate(int_exp = int_exp * sd(df_mnc$int_exp) + mean(df_mnc$int_exp))

ggplot(pd, aes(x = int_exp)) +
  geom_line(aes(y = mean, color = actorj)) +
  geom_ribbon(aes(ymin = lower95, ymax = upper95, fill = actorj), alpha = 0.25) +
  scale_color_discrete("Country") +
  scale_fill_discrete("Country") + 
  labs(x = "Export intensity", y = "Probability")
ggsave("../figure/japan96_effect_of_int_exp_on_prob_MNC_being_picked.pdf", w = 9, h = 5)

# ---- effect of log GDP on a country being chosen ----

increments <- c(0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5) # Thai gdp changes by -30%, -20%, ...
tmp <- lapply(increments, function(inc) {
  ww_hypothetical <- ww
  ww_hypothetical['Thailand', 'lgdp'] <- ww['Thailand', 'lgdp'] + log(inc)
  prob_choosing_Thai <- replicate(1000, 
    sum(f_sim_final_choice(alpha, beta, xx = xx, ww = ww_hypothetical) == "Thailand") / n_i)
  prob_choosing_Malaysia <- replicate(1000, 
    sum(f_sim_final_choice(alpha, beta, xx = xx, ww = ww_hypothetical) == "Malaysia") / n_i)
  prob_choosing_Indonesia <- replicate(1000, 
    sum(f_sim_final_choice(alpha, beta, xx = xx, ww = ww_hypothetical) == "Indonesia") / n_i)
  
  c(Thailand.mean = mean(prob_choosing_Thai), 
    Thailand = quantile(prob_choosing_Thai, probs = c(0.05, 0.95)),
    Malaysia.mean = mean(prob_choosing_Malaysia),
    Malaysia = quantile(prob_choosing_Malaysia, probs = c(0.05, 0.95)),
    Indonesia.mean = mean(prob_choosing_Indonesia),
    Indonesia = quantile(prob_choosing_Indonesia, probs = c(0.05, 0.95)))
})
pd <- as.data.frame(do.call(rbind, tmp)) %>%
  mutate(inc = increments) %>%
  reshape2::melt(id.vars = "inc") %>%
  tidyr::separate("variable", sep = "\\.", into = c("country", "var")) %>%
  spread(key = "var", value = "value")
ggplot(pd, aes(x = inc)) +
  geom_pointrange(aes(y = mean, ymin = `5%`, ymax = `95%`, col = country)) +
  scale_x_continuous(labels = scales::percent) +
  labs(x = "Thailand's hypothetical GDP", y = "Share of Japanese subsidiaries in Asia")
ggsave("../figure/japan96_effect_GDP_on_share_of_MNCs.pdf", w = 7.3, h = 4.5)

# ---- prob a country being chosen, simulate both sides ----
observed_prop <- df_mnc %>% dplyr::count(nation) %>% 
  mutate(value = n / sum(n)) %>%
  rename(choices = nation)

sim_choice <- t(replicate(1000, f_sim_final_choice(alpha, beta, xx = xx, ww = ww)))
# Coerce to data frame so that we can use dplyr::bind_rows,
# which fills missing with NA
probs_country <- apply(sim_choice, 1, 
                       function(x) as.data.frame(as.list(table(x) / n_i))) %>%
  dplyr::bind_rows() %>%
  rename(`South Korea` = "South.Korea")
f_plot_sim(probs_country) +
  theme(axis.text.x = element_text(angle = 25, vjust = 1, hjust=1)) +
  labs(x = "Country", y = "Probability") +
  geom_point(data = observed_prop, aes(x = choices, y = value, col = "Observed"))
ggsave("../figure/japan96_prob_being_chosen_by_MNC_two_sided.pdf", w = 5.5, h = 3)

# ---- How MNCs' choice change if China is more demanding? ----

simulated_location_choice <- t(replicate(1000, f_sim_final_choice(alpha, beta, xx = xx, ww = ww)))
probs_country <- apply(simulated_location_choice, 1, 
                       function(x) as.data.frame(as.list(table(x) / n_i))) %>%
  dplyr::bind_rows() %>%
  mutate(scenario = "current")

beta_china_more_demanding <- beta
beta_one <- beta[ , 'one', ]
# on avg China is now as demanding as South Korea
beta_china_more_demanding[ , 'one', "China"] <- beta_china_more_demanding[ , 'one', "China"] - 
  colMeans(beta_one)["China"] + colMeans(beta_one)["South Korea"]
simulated_location_choice2 <- t(replicate(1000, 
                   f_sim_final_choice(alpha = alpha, beta = beta_china_more_demanding, 
                                      xx = xx, ww = ww)))
probs_country2 <- apply(simulated_location_choice2, 1, 
                       function(x) as.data.frame(as.list(table(x) / n_i))) %>%
  dplyr::bind_rows() %>%
  mutate(scenario = "China demanding")
pd <- dplyr::bind_rows(probs_country, probs_country2)
molten <- reshape2::melt(pd, id.vars = "scenario", variable.name = "choices") %>%
  group_by(choices, scenario) %>%
  summarise(mean = mean(value, na.rm = TRUE), 
            lower95 = quantile(value, 0.025, na.rm = TRUE),
            upper95 = quantile(value, 0.975, na.rm = TRUE))
ggplot(molten, aes(x = choices, color = scenario)) +
  geom_pointrange(aes(y = mean, ymin = lower95, ymax = upper95)) +
  theme(axis.text.x = element_text(angle = 25, vjust = 1, hjust=1)) +
  scale_color_discrete("", 
                       labels = c("current" = "Observed situation",
                                  "China demanding" = "China becomes\nmore demanding")) +
  scale_x_discrete(labels = c("South.Korea" = "South Korea")) + 
  labs(x = "Country", y = "Probability") 
ggsave("../figure/japan96_sim_china_more_demanding.pdf", w = 7.3, h = 4.5)

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
simulated_location_choice <- t(replicate(1000, f_sim_final_choice(alpha, beta, xx = xx, ww = ww)))
simulated_profile <- apply(simulated_location_choice, 1, f_profile)

# Mean of ltemp
sim_res <- dplyr::rbind_list(lapply(simulated_profile, function(mat) mat[ , "mean(ltemp)"]))
observed_mean_ltemp <- df_mnc %>% group_by(nation) %>% 
  summarize(value = mean(ltemp)) %>%
  rename(choices = nation)

f_plot_sim(sim_res) +
  geom_point(data = observed_mean_ltemp, 
             aes(x = choices, y = value, color = "Observed")) +
  theme(axis.text.x = element_text(angle = 25, vjust = 1, hjust=1)) +
  labs(x = "Country", y = "Mean of log Employee")
ggsave("../figure/japan96_sim_mean_ltemp.pdf", w = 6, h = 3.6)

# Variance of ltemp
sim_res <- dplyr::rbind_list(lapply(simulated_profile, function(mat) mat[ , "var(ltemp)"]))
observed_var_ltemp <- df_mnc %>% group_by(nation) %>% 
  summarize(value = var(ltemp)) %>%
  rename(choices = nation)
f_plot_sim(sim_res) +
  geom_point(data = observed_var_ltemp, 
             aes(x = choices, y = value, color = "Observed")) +
  theme(axis.text.x = element_text(angle = 25, vjust = 1, hjust=1)) +
  labs(x = "Country", y = "Variance of log Employee")
ggsave("../figure/japan96_sim_var_ltemp.pdf", w = 5.5, h = 3)

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

sim_res <- t(replicate(1000, f_sim_firm_location()))
