rm(list = ls())
library("checkpoint")
checkpoint("2017-10-15", scanForPackages = FALSE)

install.packages("~/projects/match2sided/", repos = NULL, type="source")
try(detach("package:match2sided", unload = TRUE), silent = TRUE)
library("match2sided")
library("tidyverse")
library("coda")
time_start <- Sys.time()

# ---- Data ----
dat <- readRDS("../data_clean/JapanFDI_for_analysis.RDdata") %>%
  select(temp, # Total employment 
         uscptl, # Capital, in USD 
         nation_id, gdp, gdppc, democracy) %>%
  filter(uscptl > 1000) # Firms with less than 1000 USD capital is 1 pct, 
                        # likely coding error

n_i <- nrow(dat)
n_j <- length(unique(dat$nation_id))
p_i <- 3 # number of firm characteristics per firm; including the intercept
p_j <- 3  # number of country characteristics per country

# ---- Prepare choice ----

choice <- dat$nation_id

# ---- Prepare opp ----
# Opportunity sets
opp <- matrix(F, n_i, n_j)  # The opportunity matrix T=offer,F=no offer
opp[cbind(1:n_i, choice)] <- TRUE  # firms are offered the countries they are in!

# ---- Prepare ww ----
# Characteristics of country (j)

ww <- matrix(NA, n_j, p_j)
for (i in 1:n_j) {
  ind <- (1:n_i)[ dat$nation_id == i ] # index of workers whose got job i
  ww[i,1] <- unique( dat$gdp[ind] ) # gdp of nation i
  ww[i,2] <- unique( dat$gdppc[ind] ) # gdppc of nation i
  ww[i,3] <- unique( dat$democracy[ind] ) # democracy of nation i
}
# rescale (for better numerics )
ww[, 1:2] <- log(ww[, 1:2])
colnames(ww) <- c("lgdp", "lgdppc", "democracy")


# ---- Prepare xx ----

one <- rep(1, n_i)
xx <- cbind( one, dat[,1:(p_i-1)] )
# matrix n_i x p_i of firm characteristics
# including column of ones for an intercept

# rescale
xx[,2] <- xx[,2] / 10
xx[,3] <- log(xx[, 3])
xx <- as.matrix(xx)

# ---- Run MCMC ----

res <- match2sided::match2sided(iter = 10000,
                                eps_alpha = 0.05, eps_beta = 0.05, 
                                frac_beta = 0.5, frac_opp = 0.5,
                                ww = ww, xx = xx,
                                choice = choice, opp = opp)

saveRDS(res, file = paste0("../result/japan-", 
                           format(Sys.time(), "%Y-%m-%d-%H%M"), ".RData"))

# ---- Results and Diagnostics ----

WARMUP <- res$mcmc_settings$iter / 5

plot(res$lp, type='l',
     xlab = 'iteration', ylab = 'log posterior density')

alpha <- mcmc(res$alpha) %>% window(start = WARMUP)
plot(alpha)
summary(alpha)
