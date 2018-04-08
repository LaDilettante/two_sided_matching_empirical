rm(list = ls())
library("tidyverse")
library("WDI")

# ---- Global FDI ----

d_wdi_raw <- WDI(indicator = "BX.KLT.DINV.CD.WD", start = 1970, end = 2016, extra = TRUE)

d_wdi <- d_wdi_raw %>%
  filter(country %in% c("Low income", "Middle income", "High income")) %>%
  rename(fdi_inflow = BX.KLT.DINV.CD.WD) %>%
  as_tibble()

ggplot(d_wdi) +
  geom_line(aes(year, fdi_inflow, color = country)) +
  scale_y_continuous(labels = scales::unit_format("bn", scale = 1 / 1e9)) +
  scale_color_discrete("") +
  labs(x = "Year", y = "FDI inflow (USD billions)")
ggsave("../figure/global_fdi.pdf", w = 8.6, h = 4.8)

# ---- Japanese FDI ----

japFDI <- readRDS("../data_clean/japanFDI_86_99.RData")

pd <- japFDI %>%
  mutate(region = fct_collapse(region, 
                               Other = c("Oceania", "Africa", 
                                         "Middle East", "South America"))) %>%
  group_by(region, year) %>%
  summarise_at(vars(uscptl, ussale, temp), sum, na.rm = TRUE)
ggplot(pd, aes(x = year)) +
  geom_line(aes(y = temp, color = region)) +
  scale_x_continuous(breaks = 1986:1999)

japFDI %>%
  filter(region == "Asia", sic_brd == "Manufacturing") %>%
  filter(fdtn_yr > 1980) %>%
  group_by(fdtn_yr) %>%
  mutate(n = unique(sort_key))

japFDI %>%
  filter(region == "Asia", sic_brd == "Manufacturing") %>%
  filter(fdtn_yr == 1996) %>%
  summarise(n = n_distinct(sort_key))

