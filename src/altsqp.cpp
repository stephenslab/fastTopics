#include "mixsqp.h"
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// FUNCTION DECLARATIONS
// ---------------------
// TO DO.

// FUNCTION DEFINITIONS
// --------------------
// TO DO: Describe here what this function does, and describe the
// inputs and outputs.
// 
// [[Rcpp::export]]
void altsqp_update_factors_rcpp (const arma::mat& X, arma::mat& F,
				 const arma::mat& L, const arma::vec& xscol,
				 const arma::vec& ls) {
  // TO DO.
}

// TO DO: Describe here what this function does, and describe the
// inputs and outputs.
// 
// [[Rcpp::export]]
void altsqp_update_loadings_rcpp (const arma::mat& X, const arma::mat& F,
				  arma::mat& L, const arma::vec& xsrow,
				  const arma::vec& fs) {
  // TO DO.
}
