// This is included to suppress the warnings from solve() when the
// system is singular or close to singular.
#define ARMA_DONT_PRINT_ERRORS

// Disable all run-time checks, such as bounds checking. This will
// speed up code a bit.
// #define ARMA_NO_DEBUG

#include <RcppArmadillo.h>
#include "mixsqp.h"

using namespace Rcpp;
using namespace arma;

// FUNCTION DECLARATIONS
// ---------------------
// TO DO.

// FUNCTION DEFINITIONS
// --------------------
// TO DO: Describe here what this function does.
// 
// [[Rcpp::export]]
void altsqp_update_factors_rcpp (const mat& X, arma::mat& F,
				 const arma::mat& L, const arma::vec& xscol) {
  // TO DO.
}

// TO DO: Describe here what this function does.
// 
// [[Rcpp::export]]
void altsqp_update_loadings_rcpp (const mat& X, const arma::mat& F,
				  arma::mat& L, const arma::vec& xsrow) {
  // TO DO.
}
