rm(list = ls())
library("tidyverse")
library("WDI")

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
