rm(list = ls())
source("0_functions.R")

library("haven")
library("tidyverse")
library("countrycode")
library("WDI")
library("mice")

# ---- Load Japan data ----

d_raw <- readRDS("../data_clean/JapanFDI_long.RData")
d_raw_labels <- readRDS("../data_clean/JapanFDI_labels.RData")
glimpse(d_raw)

# Clean data
# Only keep countries with > 100 firms, which is about 29 countries over 139 countries
# a mix of European, NA, and SEA countries
d <- d_raw %>%
  select(id, sub_name, year, nation,
         sic_3_1, sic_2_1, jemp, temp, uscptl) %>%
  mutate(year = as.numeric(year), nation = as.character(nation),
         temp = as.numeric(temp), uscptl = as.numeric(uscptl)) %>%
  filter(year == 2001) %>% # Only keep the data in year 2001
  mutate(iso2c = countrycode(nation, origin = "country.name", destination = "iso2c")) %>%
  group_by(nation) %>% filter(n() > 100) %>% ungroup() # Keep countries with > 100 firms

# ---- Merge in country data ----

load("../data_clean/country_covariates.RData")

tmp <- d %>% distinct(iso2c, year) %>% select(iso2c, year) %>%
  left_join(d_wdi, by = c("iso2c", "year")) %>%
  left_join(d_dd, by = c("iso2c", "year"))

# Add Taiwan data
d_taiwan <- data.frame(
  iso2c = "TW", year = 2001,
  gdp = 293.68 * 1e9, gdppc = 13107.64)
d_wdi <- rbind.fill(d_wdi, d_taiwan)
d_wdi %>% filter(iso2c == "TW", year == 2001)

# Add Hong Kong data
d_hongkong <- data.frame(iso2c = "HK", year = 2001, democracy = 1)
d_dd <- rbind.fill(d_dd, d_hongkong)

d_merged <- d %>%
  left_join(d_wdi, by = c("iso2c", "year")) %>%
  left_join(d_dd, by = c("iso2c", "year")) %>%
  select(temp, uscptl, nation, gdp, gdppc, democracy)

md.pattern(d_merged %>% select(democracy, gdp, gdppc))
md.pattern(d_merged %>% select(uscptl, temp))

d_merged_final <- na.omit(d_merged)

d_nation_id <- data.frame(
  nation = sort(names(table(d_merged_final$nation))),
  nation_id = 1:length(unique(d_merged_final$nation)))

head(d_merged_final)

d_merged_final <- d_merged_final %>% inner_join(d_nation_id, by = "nation")

# ---- Create MNC dataset ----

df_mnc <- d_merged_final %>%
  mutate(luscptl = log(uscptl)) %>%
  select(temp, luscptl, country = nation, country_id = nation_id)

# ---- Create country dataset ----

df_country <- d_merged_final %>%
  mutate(lgdp = log(gdp), lgdppc = log(gdppc)) %>%
  select(country = nation, country_id = nation_id,
         lgdp, lgdppc, democracy) %>%
  distinct() %>% arrange(country_id)

# ---- Write to disk ----

write_csv(df_mnc, "../data_clean/japan_mnc.csv")
write_csv(df_country, "../data_clean/japan_country.csv")
