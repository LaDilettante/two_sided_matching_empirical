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

# Remove all countries with fewer than a certain number of firms
c_min_number_of_firms <- 1000
c_countries_to_keep <- names(table(d_raw$nation)[table(d_raw$nation) > c_min_number_of_firms])

# Clean data
d <- d_raw %>%
  select(id, sub_name, year, nation,
         sic_3_1, sic_2_1, jemp, temp, uscptl) %>%
  mutate(year = as.numeric(year), nation = as.character(nation),
         temp = as.numeric(temp), uscptl = as.numeric(uscptl)) %>%
  filter(year == 2001) %>%
  mutate(iso2c = countrycode(nation, origin = "country.name", destination = "iso2c"))

d_nation_with_more_than100firms <- d %>%
  dplyr::count(nation) %>% filter(n >= 100) %>% select(nation)

d <- d %>% inner_join(d_nation_with_more_than100firms, by = "nation")

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

# ---- Save data ----

saveRDS(d_merged_final, file = "../data_clean/JapanFDI_for_analysis.RDdata")
