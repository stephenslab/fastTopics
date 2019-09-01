#ifndef INCLUDE_MISC
#define INCLUDE_MISC

#include <RcppArmadillo.h>

// FUNCTION DECLARATIONS
// ---------------------
void   getcolnonzeros (const arma::sp_mat& A, arma::uvec& i, uint j);
void   addtorows      (arma::mat& A, const arma::vec& a);
void   scalecols      (arma::mat& A, const arma::vec& b);
void   normalizerows  (arma::mat& A);

// FUNCTION DEFINITIONS
// --------------------
// Return a or b, which ever is smaller.
inline double min (double a, double b) {
  double y;
  if (a < b)
    y = a;
  else
    y = b;
  return y;
}

#endif
