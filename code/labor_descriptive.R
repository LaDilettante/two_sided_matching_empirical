rm(list = ls())

library("tidyverse")

# ---- Data ----

# Prestige = category mean
# Autonomy = the odds of having a supervisor (multiplied by -1 so higher score = higher autonomy)
# Prestige and autonomy data not available for unemployed workers
# Prestige for unemp is set to 18, 50% of the mean prestige of the last job held
# Autonomy for unemp is set to - avg odds of having a supervisor (across the entire sample)

df_raw <- read.table("../data_clean/gss18cat.raw", header=T) %>%
  as_tibble()

df <- df_raw %>%
  mutate(occ17_label = factor(occ17, levels = 0:17, 
                              labels = c("Unemployment", 
                                         "Professionals, Self employed",
                                         "Professionals, Salaried",
                                         "Managers", "Salesmen, Other", "Proprietors",
                                         "Clerical", "Salesmen, Retail",
                                         "Craftsmen, Manufacturing",
                                         "Craftsmen, Other", "Craftsmen, Construction",
                                         "Operatives, Manufacturing", "Operatives, Other",
                                         "Service", "Laborers, Manufacturing",
                                         "Laborers, Other", "Farmers", "Farm laborers")))

unique(df %>% select(occ17_label, presmean, autmean))

hist(df$presmean)
hist(df$autmean)

lapply(df %>% select(presmean, autmean), 
       function(x) c(mean = mean(x), sd = sd(x))) %>%
  as.data.frame()

hist(df$educ)

# ---- Data viz ---

df_viz <- df %>%
  mutate(college = educ > 12)

ggplot(data = df_viz) +
  geom_histogram(aes(presmean)) +
  facet_wrap(~college, labeller = label_both) +
  labs(x = "prestige of the job")

ggplot(data = df_viz %>% filter(occ17_label %in% c("Professionals, Salaried", 
                                                   "Craftsmen, Manufacturing",
                                                   "Salesmen, Retail",
                                                   "Farm laborers"))) +
  geom_boxplot(aes(occ17_label, educ))

ggplot(data = df_viz %>% 
         filter(occ17_label %in% c("Professionals, Salaried", 
                                   "Craftsmen, Manufacturing",
                                   "Salesmen, Retail",
                                   "Farm laborers")) %>%
         select(occ17_label, presmean, autmean) %>%
         distinct()) +
  geom_col(aes(occ17_label, presmean))
