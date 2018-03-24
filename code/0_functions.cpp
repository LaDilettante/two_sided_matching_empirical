// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h> // sample.h itself #include RcppArmadillo
#include "0_functions.h"

// [[Rcpp::export]]
// Generate a vector of new offers, length n_i * number_of_new_offers
Rcpp::IntegerVector new_offer(int n_i, int n_j, int size) {
  Rcpp::IntegerVector out(n_i * size);
  
  Rcpp::IntegerVector choice = Rcpp::seq(2, n_j);
  
  for (int i = 0; i < n_i; ++i) {
    Rcpp::IntegerVector idx = Rcpp::seq(i * size, (i + 1) * size - 1);
    out[idx] = Rcpp::RcppArmadillo::sample(choice, size, false);  
  }
  
  return out;
}

// [[Rcpp::export]]
// Multivariate normal density, vector form
double dmvnrm_arma_vec(arma::vec x,  
                   arma::vec mean, arma::mat sigma, 
                   bool logd = true) { 
  
  const double log2pi = std::log(2.0 * M_PI);
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
// Multivariate normal density, matrix form
arma::vec dmvnrm_arma_mat(arma::mat x,  
                      arma::rowvec mean,  
                      arma::mat sigma, 
                      bool logd = true) { 
  const double log2pi = std::log(2.0 * M_PI);
  int n = x.n_rows;
  int xdim = x.n_cols;
  arma::vec out(n);
  arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(sigma))));
  double rootisum = arma::sum(log(rooti.diag()));
  double constants = -(static_cast<double>(xdim)/2.0) * log2pi;
  
  for (int i=0; i < n; i++) {
    arma::vec z = rooti * arma::trans( x.row(i) - mean) ;    
    out(i)      = constants - 0.5 * arma::sum(z%z) + rootisum;     
  }  
  
  if (logd == false) {
    out = exp(out);
  }
  return(out);
}

// [[Rcpp::export]]
double logsumexpC(const arma::vec& x) {
  unsigned int maxi = x.index_max();
  long double maxv = x(maxi);
  if (!(maxv > -arma::datum::inf)) {
    return -arma::datum::inf;
  }
  long double cumsum = 0.0;
  for (unsigned int i = 0; i < x.n_elem; i++) {
    if ((i != maxi) & (x(i) > -arma::datum::inf)) {
      cumsum += expl(x(i) - maxv);
    }
  }
  
  return maxv + log1p(cumsum);
}