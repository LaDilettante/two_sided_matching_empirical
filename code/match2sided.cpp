#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends("RcppArmadillo")]]

// [[Rcpp::export]]
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

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::vec dmvnrm_arma_mat(arma::mat x,  
                      arma::rowvec mean,  
                      arma::mat sigma, 
                      bool logd = false) { 
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
    dmvnrm_arma_vec(alphastar, mu, inv(Tau), true) - 
    dmvnrm_arma_vec(alpha, mu, inv(Tau), true);
  
  return logmh_alpha;
}

// logmh_beta (2410 ms for 1000 iter)
// [[Rcpp::export]]
double logmh_betaC(arma::mat beta, arma::mat betastar,
                   arma::mat xx, arma::mat opp,
                   arma::rowvec mu, arma::mat Tau) {
  arma::mat XB = xx * beta;
  arma::mat XB_star = xx * betastar;
  double lrat = accu((opp % (XB_star - XB)));
  double logmh_beta = lrat + accu(log(1 + exp(XB)) - log(1 + exp(XB_star))) +
    accu(dmvnrm_arma_mat(arma::trans(betastar), mu, inv(Tau), true)) -
    accu(dmvnrm_arma_mat(arma::trans(beta), mu, inv(Tau), true));
  return logmh_beta;
}

// logmh_opp (920 ms for 1000 iter)

// joint_logpdf (1780 ms for 1000 iter)

// f_logp_0 (660 ms for 1000 iter)

// deviation <- c(rmvnorm(1, sigma = C_ab_est)) (770 ms for 1000 iter)