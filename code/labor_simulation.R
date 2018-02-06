rm(list = ls())
library(RcppCNPy)
library(coda)
source("match2sided.R")

ww <- npyLoad("ww.npy", "numeric")
xx <- npyLoad("xx.npy", "numeric")
choice <- npyLoad("choice.npy", "integer") + 1 # Python is index 0, R is 1
opp <- npyLoad("opp.npy", "numeric")

iter <- 10000
res <- match2sided(iter = iter,
                   C_alpha = (0.2 ** 2) * diag(ncol(ww)), 
                   C_beta = (0.2 ** 2) * diag(ncol(xx)),
                   frac_opp = 0.25,
                   ww = ww, xx = xx,
                   choice = choice, opp = opp)

plot(mcmc(res$alpha))
plot(mcmc(res$alphastar))

beta_2 <- res$beta[, 2, ] # the second beta
plot(mcmc(beta_2, start = iter / 2))
