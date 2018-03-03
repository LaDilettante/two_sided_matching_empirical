#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
IntegerVector new_offer(int n_i, int n_j, int size) {
  IntegerVector out(n_i * size);
  
  IntegerVector choice = seq(2, n_j);
  
  for (int i = 0; i < n_i; ++i) {
    IntegerVector idx = seq(i * size, (i + 1) * size - 1);
    out[idx] = RcppArmadillo::sample(choice, size, false);  
  }
  
  return out;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
new_offer(10, 4, 3)
*/
