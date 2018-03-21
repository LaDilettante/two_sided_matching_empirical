#include <RcppArmadillo.h>
#include "0_functions.h"

// [[Rcpp::depends("RcppArmadillo")]]
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