rm(list = ls())
source("0_functions.R")
packs <- c("foreign", "plyr", "dplyr", "reshape2", "stringr", "tidyr",
           "lubridate", "concordance")
f_install_and_load(packs)

# ---- read data ----

d_raw <- read.spss("../raw_data/2003_entry_exit_matrix_all_data.sav",
               to.data.frame = TRUE, use.value.labels = FALSE,
               reencode = "UTF-8")
d_labels <- data.frame(name = names(d_raw),
                       label = attr(d_raw, "variable.labels"))

# ---- reshape data back to wide format ----
# The dataset is put together by merging a bunch of annual datasets
# Thus, a "minimal dataset" that preserves all the information is in wide format

d <- d_raw %>%
  select(order, sub_name,
         nation, region, province, city, zip_code,
         ind_zone, sta_pro, mode,
         fdtn_yr, fdtn_est, yeardata, yearexit, year_e,
         matches("^gm.*[0-9]{2}$"), # nationality of subsidiary general manager
         matches("^temp.*[0-9]{2}$"), # total employees
         matches("^jemp.*[0-9]{2}$"), # Japanese employees
         matches("^uscptl.*[0-9]{2}$"), # capital invested in subsidiary
         matches("^ussale.*[0-9]{2}$"), # sales of subsidiary
         matches("^prft.*[0-9]{2}$"), # profitability of subsidiary
         matches("ja1.*[0-9]{2}$"), # ownership of primary japanese partnert
         matches("^sic") # SIC code
         ) %>%
  distinct()

if (nrow(d) != length(d$order)) {
  stop("Bad: The data is not yet reduced to its minimal form.
       There are duplications that should not exist!")
} else {
  message("Good: The data is reduced to its minimal form.")
}

# ---- Convert country code
d_nation_label <- attr(d_raw$nation, "value.label")
d_nation_label <- data.frame(nation_name = names(d_nation_label),
                             nation_numeric = as.numeric(d_nation_label))

d <- d %>% inner_join(d_nation_label, by = c("nation"="nation_numeric")) %>%
  mutate(nation = nation_name) %>% select(-nation_name)

# ---- Convert industry code ----

d_sic_3_1 <- read.spss("../raw_data/2003_entry_exit_matrix_all_data.sav",
                      to.data.frame = TRUE, use.value.labels = TRUE,
                      reencode = "UTF-8") %>% select(sic_3_1)

d_industrycode <- data.frame(SIC3 = d_raw$sic_3_1,
                             SIC3_label = d_sic_3_1$sic_3_1) %>%
  distinct() %>% na.omit()
rm(d_raw, d_sic_3_1) ; gc()

tmp <- apply(d_industrycode, 1,
             function(x) concord(x[1], origin = "SIC", destination = "ISIC3"))
tmp2 <- plyr::ldply(tmp, rbind)
tmp2[, 1] <- NULL
d_industrycode <- cbind(d_industrycode, tmp2)

f_ISIC_to_intensity <- function(code) {
  library(stringr)

  hi <- "^2423|^30|^32|^33"
  medhi <- "^(?!2423)24|^29|^31|^34|^352|^359"
  medlow <- "^23|^25|^26|^27|^28|^351"
  low <- "^15|^16|^17|^19|^20|^21|^22|^36|^37"

  if (is.na(code)) {
    return(NA)
  } else if (str_detect(code, hi)) {
    return(4)
  } else if (str_detect(code, medhi)) {
    return(3)
  } else if (str_detect(code, medlow)) {
    return(2)
  } else if (str_detect(code, low)) {
    return(1)
  }
}
fv_ISIC_to_intensity <- Vectorize(f_ISIC_to_intensity)

d_industrycode2 <- d_industrycode
tmp <- lapply(d_industrycode2[, 3:18], as.character)
tmp2 <- lapply(lapply(d_industrycode2[, 3:18], as.character), fv_ISIC_to_intensity)
d_industrycode2[, 3:18] <- data.frame(do.call(cbind, tmp2))
is.na(d_industrycode2) <- d_industrycode2 == "NULL"
d_industrycode2[, 3:18] <- lapply(d_industrycode2[, 3:18], as.numeric)
d_industrycode2$intensity_avg <- rowMeans(d_industrycode2[ , 3:18], na.rm = TRUE)
d_industrycode2$intensity_avg[is.nan(d_industrycode2$intensity_avg)] <- NA

# Merge industry intensity into main data

attr(d$sic_3_1, "value.labels") <- NULL
d <- left_join(d, d_industrycode2[, c("SIC3", "SIC3_label", "intensity_avg")],
               by = c("sic_3_1"="SIC3"))

# ---- Save wide data ----

saveRDS(d, file = "../clean_data/JapanFDI_wide.RData")
saveRDS(d_labels, file = "../clean_data/JapanFDI_labels.RData")
