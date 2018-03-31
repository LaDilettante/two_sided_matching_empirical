#include <RcppArmadillo.h>
#include "0_functions.h"
// [[Rcpp::depends("RcppArmadillo")]]

// [[Rcpp::export]]
// This version is faster, but not safe for numerical overflow
// due to exp
double logmh_alphaC_old(arma::vec alpha, arma::vec alphastar,
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

// [[Rcpp::export]]
double logmh_alphaC(arma::vec alpha, arma::vec alphastar,
                   arma::mat ww, arma::mat oppT,
                   arma::vec wa,
                   arma::vec mu, arma::mat Tau) {
  // Calculate likelihood ratio
  int n_i = oppT.n_cols;
  arma::vec walpha = ww * alpha;
  arma::vec pA_den(n_i);
  for (arma::uword i = 0; i < n_i; i++) {
    pA_den(i) = logsumexpC(walpha.elem(arma::find(oppT.col(i))));
  }
  
  arma::vec walphastar = ww * alphastar;
  arma::vec pA_denstar(n_i);
  for (arma::uword i = 0; i < n_i; i++) {
    pA_denstar(i) = logsumexpC(walphastar.elem(arma::find(oppT.col(i))));
  }
  
  double logmh_alpha = sum(wa % (alphastar - alpha)) + 
    sum(pA_den - pA_denstar) +
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
  double p2 = sum(Rcpp::log1p(Rcpp::wrap(exp(XB))) - Rcpp::log1p(Rcpp::wrap(exp(XB_star))));
  double p3 = accu(dmvnrm_arma_mat(arma::trans(betastar), mu, inv(Tau), true)) -
    accu(dmvnrm_arma_mat(arma::trans(beta), mu, inv(Tau), true));
  
  double logmh_beta = lrat + p2 + p3;
  return logmh_beta;
}

// logmh_opp (920 ms for 1000 iter)


// [[Rcpp::export]]
double f_logp_AC(arma::mat opp, arma::vec alpha, 
                 arma::mat ww, arma::mat w_chosen) {
  return(accu(w_chosen * alpha) - accu(log(opp * exp(ww * alpha))));
}

// f_logp_0 (660 ms for 1000 iter)
// [[Rcpp::export]]
double f_logp_OC(arma::mat opp, arma::mat beta, arma::mat xx) {
  arma::mat XB = xx * beta;
  return(accu(opp % XB) - accu(log(1 + exp(XB))));
}

// joint_logpdf (1780 ms for 1000 iter)

// deviation <- c(rmvnorm(1, sigma = C_ab_est)) (770 ms for 1000 iter)