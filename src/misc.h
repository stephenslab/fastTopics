#ifndef INCLUDE_MISC
#define INCLUDE_MISC

#include <RcppArmadillo.h>

// FUNCTION DECLARATIONS
// ---------------------
double min           (double a, double b);
void   getnonzeroindicesincol (const arma::sp_mat& A, arma::uvec& i, uint j);
void   addtorows     (arma::mat& A, const arma::vec& a);
void   scalecols     (arma::mat& A, const arma::vec& b);
void   normalizerows (arma::mat& A);
#endif
