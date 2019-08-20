#ifndef INCLUDE_MISC
#define INCLUDE_MISC

#include <RcppArmadillo.h>

// FUNCTION DECLARATIONS
// ---------------------
double min          (double a, double b);
void   copycolelems (const arma::mat& X, const arma::uvec& i, uint j,
		     arma::vec& y);
void   copyrowelems (const arma::mat& X, uint i, const arma::uvec& j,
		     arma::vec& y);
void   scalecols    (arma::mat& A, const arma::vec& b);
  
#endif
