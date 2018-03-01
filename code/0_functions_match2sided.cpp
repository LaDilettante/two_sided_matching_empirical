#include <RcppArmadillo.h>
using namespace Rcpp;

const double log2pi = std::log(2.0 * M_PI);

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
double dmvnrm_arma(arma::vec x,  
                   arma::vec mean, arma::mat sigma, 
                   bool logd = true) { 
  double out;
  int xdim = x.n_elem;
  arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(sigma))));
  double rootisum = arma::sum(log(rooti.diag()));
  double constants = -(static_cast<double>(xdim)/2.0) * log2pi;
  
  arma::vec z = rooti * (x - mean) ;    
  out = constants - 0.5 * arma::sum(z%z) + rootisum;     
  
  if (logd == false) {
    out = exp(out);
  }
  return(out);
}

// [[Rcpp::export]]
double logmh_alphaC(arma::vec alpha, arma::vec alphastar,
                   arma::mat ww, arma::mat opp,
                   arma::vec wa,
                   arma::vec mu, arma::mat Tau) {
  // Calculate likelihood ratio
  arma::vec exp_WA = exp(ww * alpha);
  arma::vec pA_den = opp * exp_WA;
  arma::vec exp_WA_star = exp(ww * alphastar);
  arma::vec pA_denstar = opp * exp_WA_star;
  
  double logmh_alpha = sum(wa % (alphastar - alpha)) + 
    sum(log(pA_den) - log(pA_denstar)) +
    dmvnrm_arma(alphastar, mu, inv(Tau), true) - 
    dmvnrm_arma(alpha, mu, inv(Tau), true);
  
  return logmh_alpha;
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
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

# Test
*/
