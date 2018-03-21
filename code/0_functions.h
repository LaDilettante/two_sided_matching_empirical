#ifndef GUARD_functions_h
#define GUARD_functions_h

#include <RcppArmadillo.h>

Rcpp::IntegerVector new_offer(int n_i, int n_j, int size);

double dmvnrm_arma_vec(arma::vec x, arma::vec mean, arma::mat sigma, 
                       bool logd);

arma::vec dmvnrm_arma_mat(arma::mat x,  
                         arma::rowvec mean,  
                         arma::mat sigma, 
                         bool logd);
  
#endif