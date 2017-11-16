rm(list = ls())
source("0_functions.R")

library("tidyverse")
library("ggplot2")
library("countrycode")

# ---- Load Japan data ----

d_raw <- readRDS("../data_clean/JapanFDI_long.RData")
d_raw_labels <- readRDS("../data_clean/JapanFDI_labels.RData")
glimpse(d_raw)

table(d_raw$nation)
