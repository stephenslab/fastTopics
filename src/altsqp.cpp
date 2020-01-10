#include <RcppParallel.h>
#include "misc.h"
#include "poismix.h"

using namespace arma;

// INLINE FUNCTION DEFINITIONS
// ---------------------------
inline vec altsqp_update_factor (const mat& X, const mat& F, const mat& L1,
				 const vec& u, mat& Z, mat& H, uint j,
				 uint numiter,
				 const mixsqp_control_params& control) {
  vec f = F.col(j);
  poismixsqp(L1,u,X.col(j),f,Z,H,numiter,control);
  return f;
}

// CLASS DEFINITIONS
// -----------------
// TO DO.

// FUNCTION DEFINITIONS
// --------------------
// TO DO: Explain here what this function does, and how to use it.
// 
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat altsqp_update_factors_rcpp (const arma::mat& X, const arma::mat& F,
				      const arma::mat& L, double numiter,
				      const Rcpp::List control) {
  mixsqp_control_params ctrl = get_mixsqp_control_params(control);
  uint m    = X.n_cols;
  vec  u    = sum(L,0);
  mat  L1   = L;
  mat  Z    = L;
  mat  Fnew = F;
  mat  H(m,m);
  normalizecols(L1);
  for (uint j = 0; j < m; j++)
    Fnew.col(j) = altsqp_update_factor(X,F,L1,u,Z,H,j,numiter,ctrl);
  return Fnew;
}

// This does the same thing as altsqp_update_factors_rcpp, except that
// X is a sparse matrix. See altsqp_update_factors_rcpp for details.
//
// [[Rcpp::export]]
arma::mat altsqp_update_factors_sparse_rcpp (const arma::sp_mat& X,
					     const arma::mat& F,
					     const arma::mat& L,
					     double numiter) {
  mat Fnew = F;
  return Fnew;
}
