#ifndef INCLUDE_MISC
#define INCLUDE_MISC

#include <RcppArmadillo.h>

#define maximum(a,b) ((a) > (b) ? (a) : (b))
#define minimum(a,b) ((a) < (b) ? (a) : (b))

// FUNCTION DECLARATIONS
// ---------------------
void   getcolnonzeros      (const arma::sp_mat& A, arma::uvec& i, uint j);
void   scalecols           (arma::mat& A, const arma::vec& b);
void   normalizerows       (arma::mat& A);
void   normalizecols       (arma::mat& A);
void   normalizerowsbymax  (arma::mat& A);
double dot_sparse_b        (const arma::vec& a, const arma::vec& b,
			    const arma::uvec& i);
double dot_square_sparse_b (const arma::vec& a, const arma::vec& b,
			    const arma::uvec& i);

#endif
