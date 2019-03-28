// This is included to suppress the warnings from solve() when the
// system is singular or close to singular.
#define ARMA_DONT_PRINT_ERRORS

#include <RcppArmadillo.h>

// This depends statement is needed to tell R where to find the
// additional header files.
//
// [[Rcpp::depends(RcppArmadillo)]]
//

using namespace Rcpp;

// TO DO: Explain here what this function does, and how to use it.
// 
// [[Rcpp::export]]
//
double activeset_rcpp (const arma::mat& H, const arma::vec& g,
		       const arma::vec& x0) {
  return g.min();
}
