im(list = ls())
source("0_functions.R")

# library("checkpoint")
# checkpoint("2017-10-15")

library("tidyverse")
library("stringr")
library("lubridate")
library("countrycode")
library("WDI")

# ---- Read data ----

# NATION has no value label for NATION == 200, but it's ok since there's only 1 firm
# Al-Ghalilah Engineering & Construction L
# Duplicated levels in factor PAR_AGE2: is coerced to factor because it has a value label for 0, which means that sub is founded earlier than parent
# Undeclared level(s) 0.579999999999998, 0.670000000000002
# Undeclared level(s) 5, 18, 43, 69, 71 in SIC_2 (some values that don't exist in SIC_2)
d_raw <- foreign::read.spss("../data_raw/1999_97_94_92_86 all variables.sav",
                            to.data.frame = TRUE, use.value.labels = TRUE,
                            reencode = "UTF-8")
d_labels <- data.frame(name = names(d_raw),
                       label = attr(d_raw, "variable.labels"))

d_wide <- d_raw %>%
  rename_all(tolower) %>%
  select(sort_key, merge, sub_name,  
         nation, region,
         firmname, p.1, par_list, par_ind1,
         mode,
         fdtn_dec, fdtn_yr, sub_age, par_age1, par_herf, 
         sic_2, sic_nrw, sic_brd, sur_ext, yeardata, yearexit,
         keiretsu, 
         starts_with("exp_"),
         sale_net, sale_ln,
         employee, emp_ln, size_re1, size_re2,
         starts_with("int_"), income_o, income_n,
         matches("^gm.*[0-9]{2}$"), # nationality of subsidiary general manager
         matches("^temp.*[0-9]{2}$"), # total employees
         matches("^jemp.*[0-9]{2}$"), # Japanese employees
         matches("^uscptl.*[0-9]{2}$"), # capital invested in subsidiary
         matches("^ussale.*[0-9]{2}$"), # sales of subsidiary
         matches("^prft.*[0-9]{2}$"), # profitability of subsidiary
         matches("ja1.*[0-9]{2}$") # ownership of primary japanese partnert
  ) %>%
  distinct() %>% as_tibble()

# ---- Clean nation label ----

d_wide <- d_wide %>%
  filter(!nation %in% c("Virgin Islands", "Saipan Island", 
                        "Netherlands Antilles", "200", "Dutch Antilles")) %>%
  mutate(nation = countrycode(nation, "country.name", "country.name",
                              custom_match = c("Czech" = "Czech Republic"))) %>%
  mutate(iso3c = countrycode(nation, "country.name", "iso3c"))

# ---- Reshape from wide to long ----

# Time-invariant variables, i.e. those not ending in years
c_id_vars <- names(d_wide)[!stringr::str_detect(names(d_wide), "86$|92$|94$|97$|99$")]

d_molten <- reshape2::melt(d_wide, id.vars = c_id_vars)

# Split the variable into variable and year
d_molten <- separate(d_molten, col = "variable", into = c("variable2", "year"),
                     sep = "(?<=ja1)_|(?<=\\w)(?=86|92|94|97|99)",
                     remove = F, extra = "merge")

# Check if the split works
unique(d_molten[, c("variable", "variable2", "year")])
d_molten <- d_molten %>% select(-variable)

# Reformat the year 2 digits to 4 digits
d_molten$year <- parse_date_time(d_molten$year, orders = "y")
d_molten$year <- format(d_molten$year, "%Y")

# Cast the data & and select the variables

# Note that year is the last variable on the left hand side,
# and thus is the variable that varies the fastest
fm_cast <- as.formula(str_c(str_c(c_id_vars, collapse = " + "), " + year ~ variable2"))
d_cast <- reshape2::dcast(d_molten, formula = fm_cast) %>%
  select(sort_key, sub_name, p.1, firmname, nation, year, everything()) %>%
  mutate_at(vars(year, ja1, temp, prft, temp, uscptl, ussale), funs(as.numeric)) %>%
  mutate(ltemp = log(temp), luscptl = log(uscptl), lussale = log(ussale)) %>%
  arrange(sub_name, firmname, year) %>%
  as_tibble()

saveRDS(d_cast, file = "../data_clean/japanFDI_86_99.RData")

# ---- Clean country data ----

d_pwt <- haven::read_dta("../data_raw/pwt90.dta")
d_pwt_label <- tibble(var = names(d_pwt),
                      label = sapply(d_pwt, attr, "label"))

d_country <- d_pwt %>%
  select(country, countrycode, year, # ISO-3
         pop, # population, millions
         emp, # labor force
         hc, # human capital index
         rgdpo, # output real GDP, chained PPP, 2011 USD mil
         labsh, # share of labor compensation in GDP at current national prices
         csh_x, csh_m) %>% # export / import as share of real GDP, current PPP
  mutate(lpop = log(pop),
         lgdp = log(rgdpo),
         lgdppc = log(rgdpo / pop),
         labor_cost = labsh * rgdpo / emp,
         pct_trade_share_gdp = csh_x + csh_m) %>%
  group_by(country) %>% arrange(country, year) %>%
  mutate(emp_growth = (emp - dplyr::lag(emp, 1)) / dplyr::lag(emp, 1)) %>%
  ungroup() %>%
  # Cut down the samples
  filter(year >= min(d_cast$year), year <= max(d_cast$year)) %>%
  filter(countrycode %in% unique(d_cast$iso3c)) %>%
  # Cut down the variables
  select(country, iso3c = countrycode, year,
         lpop, lgdp, lgdppc,
         hc, emp, emp_growth, labor_cost, pct_trade_share_gdp)

# Labor force data from WDI
d_labor_raw <- WDI::WDI(country = countrycode(unique(d_country$iso3c), "iso3c", "iso2c"),
                        indicator = c("SL.TLF.TOTL.IN"),
                        start = min(d_cast$year) - 1, end = max(d_cast$year))
d_labor <- d_labor_raw %>%
  select(iso2c, year, labor_force = SL.TLF.TOTL.IN) %>%
  mutate(country = countrycode(iso2c, "iso2c", "country.name")) %>%
  arrange(country, year) %>%
  group_by(country) %>%
  mutate(labor_force_growth = ((labor_force - dplyr::lag(labor_force, 1)) / dplyr::lag(labor_force, 1)))

# Labor force data from ILO
d_ilo_raw <- read_csv("../data_raw/ILO_wage.csv")
d_ilo <- d_ilo_raw %>%
  filter(classif1.label == "Activity: Manufacturing",
         sex.label == "Sex: Total") %>%
  select(country = ref_area.label,
         year = time, value = obs_value,
         sex, indicator, indicator.label) %>%
  mutate(country = countrycode(country, "country.name", "country.name"),
         iso3c = countrycode(country, "country.name", "iso3c")) %>%
  mutate(year = as.numeric(year)) %>%
  # Cut down the samples
  filter(year >= min(d_cast$year), year <= max(d_cast$year)) %>%
  filter(iso3c %in% unique(d_cast$iso3c)) %>%
  # Cut down the variables
  select(iso3c, year, avg_wage = value)

saveRDS(d_country, file = "../data_clean/country_86_99.RData")
