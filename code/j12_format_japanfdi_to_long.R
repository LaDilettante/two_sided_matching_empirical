rm(list = ls())
source("0_functions.R")
packs <- c("dplyr", "reshape2", "stringr", "tidyr",
           "lubridate")
f_install_and_load(packs)

# ---- Load the clean wide data ----

d <- readRDS("../data_clean/JapanFDI_wide.RData")

# ---- reshape the clean wide format into long format ----

# Time-invariant variables
c_id_vars <- c("order", "sub_name",
               "nation", "region", "province", "city", "zip_code", "ind_zone",
               "sta_pro", "mode", "fdtn_yr", "fdtn_est", "yeardata", "yearexit",
               "year_e", "sic_3_1", "sic_3_2", "sic_3_3",
               "sic_3_4", "sic_3_5", "sic_2_1", "sic_2_2", "sic_2_3", "sic_2_4",
               "sic_2_5", "sic_brd1", "sic_brd2", "sic_brd3", "sic_brd4",
               "sic_brd5", "sic_nrw1", "sic_nrw2", "sic_nrw3", "sic_nrw4",
               "sic_nrw5", "SIC3_label", "intensity_avg")

d_molten <- melt(d, id.vars = c_id_vars)

# Split the variable into variable and year
d_molten <- separate(d_molten, col = "variable", into = c("variable2", "year"),
                     sep = "(?<=ja1)_|(?<=\\w)(?=86|92|94|97|99|01|03)",
                     remove = F, extra = "merge")

# Check if the split works
unique(d_molten[, c("variable", "variable2", "year")])
d_molten <- d_molten %>% select(-variable)

# Reformat the year 2 digits to 4 digits
d_molten$year <- parse_date_time(d_molten$year, orders = "y")
d_molten$year <- format(d_molten$year, "%Y")

# Cast the data

# Note that year is the last variable on the left hand side,
# and thus is the variable that varies the fastest
fm_cast <- as.formula(str_c(str_c(c_id_vars, collapse = " + "), " + year ~ variable2"))
d_cast <- dcast(d_molten, formula = fm_cast)
d_cast <- d_cast %>%
  select(order, sub_name, year, everything()) %>%
  arrange(order, year)

# ---- Final cleaning ----
rm(d_molten) ; gc()

d_cast <- tbl_df(d_cast) %>%
  select(id = order, everything())

# ---- Save the data ----
saveRDS(d_cast, file = "../data_clean/JapanFDI_long.RData")
