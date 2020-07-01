#ifndef INCLUDE_MISC
#define INCLUDE_MISC

#include <RcppArmadillo.h>

#define maximum(a,b) ((a) > (b) ? (a) : (b))
#define minimum(a,b) ((a) < (b) ? (a) : (b))

// FUNCTION DECLARATIONS
// ---------------------
void getcolnonzeros     (const arma::sp_mat& A, arma::uvec& i, unsigned int j);
void getcolelems        (const arma::mat& A, const arma::uvec& i, 
			 unsigned int j, arma::vec& a);
void scalecols          (arma::mat& A, const arma::vec& b);
void normalizerows      (arma::mat& A);
void normalizecols      (arma::mat& A);
void normalizerowsbymax (arma::mat& A);

#endif
