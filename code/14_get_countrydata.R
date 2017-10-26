rm(list = ls())
source("0_functions.R")
packs <- c("haven", "reshape2", "plyr", "dplyr", "WDI",
           "countrycode", "zoo", "stringr")
f_install_and_load(packs)

# Note: Merge everything using iso2c

# ---- Constants ----
c_startyear = 1995
c_endyear = 2010

# ---- Get WDI data ----

d_wdi_raw <- WDI(
  country = "all",
  indicator = c("NY.GDP.MKTP.KD", "NY.GDP.PCAP.KD", "SL.TLF.TERT.ZS",
                "SE.SEC.NENR", "SL.TLF.SECO.ZS"),
  start = c_startyear, end = c_endyear, extra = TRUE)

d_wdi <- d_wdi_raw %>%
  filter(region != "Aggregates") %>%
  select(iso2c, country, year,
         gdp = NY.GDP.MKTP.KD,
         gdppc = NY.GDP.PCAP.KD,
         enroll_sec_pct = SE.SEC.NENR,
         labor_sec_pct = SL.TLF.SECO.ZS) %>%
  arrange(country, year)

# Supplement with Singapore data

d_sing <- read.csv("../raw_data/singapore/net-enrolment-ratio-primary-and-secondary-education.csv") %>%
  filter(gender == "both sexes", level_of_education == "secondary education") %>%
  mutate(country = "Singapore") %>%
  select(country, year, net_enrolment_ratio)

d_wdi <- left_join(d_wdi, d_sing, by = c("country", "year"))
d_wdi <- d_wdi %>%
  mutate(enroll_sec_pct = ifelse(is.na(enroll_sec_pct), net_enrolment_ratio, enroll_sec_pct)) %>%
  select(-net_enrolment_ratio)

# Supplement with education data

d_avgschooling <- read.csv("../raw_data/mean_year_of_schooling.csv",
                           na.strings = "..")
d_avgschooling$X2003 <- (d_avgschooling$X2001 + d_avgschooling$X2005) / 2
d_avgschooling <- melt(d_avgschooling, id.vars = "Country", variable.name = "year") %>%
  mutate(year = as.numeric(str_sub(year, start = 2)),
         avg_schooling_years = value,
         iso2c = countrycode(Country, origin = "country.name",
                             destination = "iso2c")) %>%
  select(iso2c, year, avg_schooling_years)

d_wdi <- left_join(d_wdi, d_avgschooling, by = c("iso2c", "year"))


# ---- Get DPI data ----

d_dpi_raw <- psData::DpiGet()
d_dpi <- d_dpi_raw %>%
  select(iso2c, countryname, year, system, gov1age, execage)

# ---- Clean dd data ----

d_dd_raw <- read_dta("../raw_data/GWF Autocratic Regimes 1.2/GWF_AllPoliticalRegimes.dta")

d_dd <- d_dd_raw %>%
  filter(year >= c_startyear, year <= c_endyear) %>%
  mutate(iso2c = countrycode(cowcode, origin="cown", destination="iso2c", warn=TRUE)) %>%
  mutate(democracy = ifelse(gwf_nonautocracy == "democracy", 1, 0)) %>%
  select(iso2c, country = gwf_country, year, democracy, spell = gwf_spell, cowcode)

# ---- Save data ----

save(d_dd, d_wdi, d_dpi, file = "../clean_data/country_covariates.RData")

