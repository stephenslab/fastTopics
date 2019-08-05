#include <RcppArmadillo.h>

// This is needed to tell R where to find the additional header files.
//
// [[Rcpp::depends(RcppArmadillo)]]
//

using namespace Rcpp;
using namespace arma;

// FUNCTION DECLARATIONS
// ---------------------

// FUNCTION DEFINITIONS
// --------------------
// [[Rcpp::export]]
double altsqp_main_loop_rcpp (const arma::mat& X, arma::mat& F, arma::mat& Fn,
			      arma::mat& Fy, arma::mat& Fbest, arma::mat& L,
			      arma::mat& Ln, arma::mat& Ly, arma::mat& Lbest,
			      double f, double fbest, bool verbose) {

  // xsrow,xscol,beta,betamax,numiter,control,progress,verbose)

  if (verbose)
    Rprintf("Running altsqp_main_loop_cpp\n");

  return fbest;
}
