rm(list = ls())

install.packages("~/projects/match2sided/", repos = NULL, type="source")
detach("package:match2sided", unload = TRUE)
library("match2sided")
library("tidyverse")
library("coda")
time_start <- Sys.time()

# ---- Data ----
dat <- read.table("../data_clean/gss18cat.raw", header=T)
choice <- dat$occ17 + 1 # so 1=unemployment
n_i <-  length(choice) # number_of_employees, Job acceptances (elements in 1:n_j)

n_j <- length(unique(choice))    # number of employers, includes unemployment
p_i <- 4 # number of worker characteristics per worker; including the intercept
# in this example it measures quality of employee
p_j <- 2  # number of job characteristics per job; in this example, job quality

# Populate the opportunity set
opp <- matrix(F, n_i, n_j)  # The opportunity matrix T=offer,F=no offer
opp[cbind(1:n_i, choice)] <- T  # people are offered jobs they have!
opp[, 1] <- T                     # Unemployment always offered

# Populate the W_j matrix
ww <- matrix(NA, n_j, p_j)
for (i in 1:n_j) {
  ww[i, 1] <- unique(dat$presmean[choice == i])
  ww[i, 2] <- unique(dat$autmean[choice == i])
}
ww <- scale(ww) # rescale (for better numerics )
colnames(ww) <- c("prestige", "autonomy")
rownames(ww) <- c("Unemployment", "Professionals, Self employed",
                  "Professionals, Salaried",
                  "Managers", "Salesmen, Other", "Proprietors",
                  "Clerical", "Salesmen, Retail",
                  "Craftsmen, Manufacturing",
                  "Craftsmen, Other", "Craftsmen, Construction",
                  "Operatives, Manufacturing", "Operatives, Other",
                  "Service", "Laborers, Manufacturing",
                  "Laborers, Other", "Farmers", "Farm laborers")

# Populate the X_i matrix
# matrix nworkers x nx of worker characteristics
# including column of ones for an intercept
one <- rep(1, n_i)
xx <- scale(dat[, c("educ", "age", "nonwhite")])
xx <- cbind(one, xx)

# ---- Run MCMC ----

res <- match2sided::match2sided(iter = 10000,
                  eps_alpha = 0.05, eps_beta = 0.05, 
                  frac_beta = 0.5, frac_opp = 0.5,
                  ww = ww, xx = xx,
                  choice = choice, opp = opp)

saveRDS(res, file = paste0("../result/labor-", 
                           format(Sys.time(), "%Y-%m-%d-%H%M"), ".RData"))

# ---- Results and Diagnostics ----

WARMUP <- res$mcmc_settings$iter / 5

plot(res$lp, type='l',
     xlab = 'iteration', ylab = 'log posterior density')

alpha <- mcmc(res$alpha) %>% window(start = WARMUP)
plot(alpha)
summary(alpha)

beta_educ <- mcmc(res$beta[, 'educ', ]) %>% window(start = WARMUP)
plot(beta_educ[, c('Professionals, Salaried', 'Farm laborers')])
plot(beta_educ)

beta_age <- mcmc(res$beta[, 'age', ]) %>% window(start = WARMUP)
plot(beta_age[, c('Professionals, Salaried', 'Farm laborers')])
plot(beta_age)



print(res$acceptance_rate)


boxplot(educ ~ occ17, data = dat %>% filter(occ17 %in% c(2, 3, 10, 17)),
        names = c("professional", "managers",
                  "construction worker", "farm worker"),
        ylab = "years of education",
        main = "sample statistics of education")


par(mfrow = c(1, 2))
boxplot(educ ~ occ17, data = dat %>% filter(occ17 %in% c(2, 17)),
        names = c("professional", "farm worker"),
        ylab = "years of education",
        main = "sample statistics of education")

plot(density(bsave[500:mcmc$nsave, 4 * (2 - 1) + 1]),
     xlim = c(-0.5, 0.5),
     main = "estimated firms' preference for education")
lines(density(bsave[500:mcmc$nsave, 4 * (17 - 1) + 1]),
      lty = 2)
legend(-0.45, 5, c("profressional", "farm worker"), lty = c(1, 2))
par(mfrow = c(1, 1))

# ---- Data viz ---

viz_dat <- dat
viz_dat <- viz_dat %>%
  mutate(college = educ > 12)

ggplot(data = viz_dat) +
  geom_density(aes(presmean)) +
  facet_wrap(~college, labeller = label_both) +
  labs(x = "prestige of the job")