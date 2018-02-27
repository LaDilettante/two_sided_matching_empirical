rm(list = ls())
library("tidyverse")
library("forcats")

# ---- Data ----
df_raw <- read.table("../data_clean/gss18cat.raw", header=T)
df <- df_raw %>%
  as_tibble() %>%
  dplyr::rename(occ17_num = occ17) %>%
  mutate(occ17 = factor(occ17_num))

# Labels come from email of Newton / Logan
levels(df$occ17) <- c("Unemployment", "Professionals, Self employed",
                       "Professionals, Salaried",
                       "Managers", "Salesmen, Other", "Proprietors",
                       "Clerical", "Salesmen, Retail",
                       "Craftsmen, Manufacturing",
                       "Craftsmen, Other", "Craftsmen, Construction",
                       "Operatives, Manufacturing", "Operatives, Other",
                       "Service", "Laborers, Manufacturing",
                       "Laborers, Other", "Farmers", "Farm laborers")

# Manually match up based on Newton email and compare with AJS paper
# Newton email has some mistakes
df$occ5 <- fct_collapse(df$occ17,
  "Professional" = c("Professionals, Self employed", "Professionals, Salaried"),
  "Managerial" = c("Managers", "Salesmen, Other", "Proprietors", "Farmers"),
  "Sales, Clerical, Services" = c("Clerical", "Salesmen, Retail", "Service"),
  "MFG Blue Collar" = c("Craftsmen, Manufacturing", "Operatives, Manufacturing", 
                        "Laborers, Manufacturing"),
  "Other Blue Collar" = c("Craftsmen, Other", "Craftsmen, Construction",
                          "Operatives, Other", "Laborers, Other", "Farm laborers"))

# ---- Create employee level dataset ----

df_employee <- df %>%
  select(educ, age, nonwhite, occ17_num, occ17, occ5)
write_csv(df_employee, "../data_clean/labor_employee.csv")

# ---- Create the reduced dataset as in AJS ----

df_occ17 <- df %>%
  select(occ17, presmean, autmean) %>%
  distinct() %>%
  mutate(supervisor = (- autmean) / (1 - autmean))
write_csv(df_occ17, "../data_clean/labor_employer_occ17.csv")

df_occ5 <- df %>% 
  select(occ17, occ5, presmean, autmean) %>%
  group_by(occ17, occ5) %>% # Collapse individual level to occ17 level
  summarise(n = n(), pres = unique(presmean), aut = unique(autmean)) %>%
  mutate(supervisor = (- aut) / (1 - aut)) %>% # Create P(supervisor) from odd
  group_by(occ5) %>% # Collapse occ17 level to occ5 level
  summarise(n_obs = sum(n),
            pres = weighted.mean(pres, n),
            aut = weighted.mean(aut, n),
            supervisor = weighted.mean(supervisor, n)) %>%
  select(-n_obs)
write_csv(df_occ5, "../data_clean/labor_employer_occ5.csv")