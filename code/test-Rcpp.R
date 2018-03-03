rm(list = ls())

context("Metropolis-Hasting acceptance ratio Rcpp")
source("0_functions.R")
source("match2sided.R")
Rcpp::sourceCpp("match2sided.cpp")
Rcpp::sourceCpp("0_functions.cpp")

# ---- Test helper functions ----

test_that("dmvnrm_arma calculates the right MVN density", {
  dim <- sample(2:10, size = 1) # a random dimension
  mu <- rnorm(dim)
  Sigma <- diag(abs(rnorm(dim)))
  x <- rmvnorm(1, mu, Sigma)
  expected <- mvtnorm::dmvnorm(x, mu, Sigma, log = TRUE)
  observed <- dmvnrm_arma(x, mu, Sigma, TRUE)
  expect_equal(observed, expected)
})

test_that("new_offer generates the right offer", {
  n_i <- sample(2000:3000, size = 1)
  n_j <- sample(10:30, size = 1)
  size <- sample(1:n_j, size = 1)
  set.seed(1)
  expected <- c(replicate(n_i, sample(2:n_j, size = size, replace = FALSE)))
  set.seed(1)
  observed <- new_offer(n_i, n_j, size)
  expect_equal(expected, observed)
})

# ---- Simulate data ----
n_i <- sample(seq(10, 20), 1) # multinomial logit side
p_i <- sample(seq(2, 3), 1)
n_j <- sample(seq(5, 10), 1) # logit side
p_j <- sample(seq(2, 3), 1)

choice <- sample(seq(1, n_j), size = n_i, replace = TRUE) # Each i picks a j
opp <- matrix(FALSE, n_i, n_j)
opp[cbind(1:n_i, choice)] <- TRUE # the choice is always in the opportunity set
# randomly flip 1/10 to TRUE
opp[sample(1:(n_i * n_j), size = n_i * n_j / 10, replace = TRUE)] <- TRUE

alpha <- rnorm(p_j)
beta <- matrix(rnorm(p_i * n_j), p_i, n_j)
mu_beta <- rnorm(p_i)
Tau_beta <- solve(diag(abs(rnorm(p_i))))
prior <- list(alpha = list(mu = rnorm(p_j), Tau = solve(diag(abs(rnorm(p_j))))),
              mu_beta = list(mu = rnorm(p_i), Tau = solve(diag(abs(rnorm(p_i))))),
              Tau_beta = list(nu = p_i + 2, Sinv = solve(diag(abs(rnorm(p_i))))))

ww <- matrix(rnorm(n_j * p_j), n_j, p_j)
tmp <- as.matrix(ww[choice, ])
wa <- apply(tmp, 2, sum)

xx <- matrix(rnorm(n_i * p_i), n_i, p_i) # Characteristics of the multinomial side
xx[, 1] <- 1 # The first column is the intercept

# ---- Tests ----

test_that("logmh_alpha is correct", {
  eps <- rnorm(1)
  deviation <- eps * runif(length(alpha), min=-1, max=1) # Symmetric proposal
  alphastar <- alpha + deviation
  
  expected <- logmh_alpha(alpha, alphastar, ww, opp, wa, prior)
  observed <- logmh_alphaC(alpha, alphastar, ww, opp, wa, prior$alpha$mu, prior$alpha$Tau)
  expect_equal(expected, observed)
})

# ---- Benchmark ----

eps <- rnorm(1)
deviation <- eps * runif(length(alpha), min=-1, max=1) # Symmetric proposal
alphastar <- alpha + deviation
microbenchmark::microbenchmark(
  logmh_alpha(alpha, alphastar, ww, opp, wa, prior),
  logmh_alphaC(alpha, alphastar, ww, opp, wa, prior$alpha$mu, prior$alpha$Tau)
)


n_i <- 2000
n_j <- 20
size <- 10
microbenchmark::microbenchmark(
  replicate(n_i, sample(2:n_j, size = size, replace = FALSE)),
  new_offer(n_i, n_j, size)
)


