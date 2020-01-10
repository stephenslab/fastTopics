#include <RcppParallel.h>
#include "misc.h"
#include "poismix.h"

using namespace arma;

// INLINE FUNCTION DEFINITIONS
// ---------------------------
// Perform one or more alt-SQP updates for a single column of the k x m
// factors matrix, F, when X is a dense matrix.
inline vec altsqp_update_factor (const mat& X, const mat& F, const mat& L1,
				 const vec& u, mat& Z, mat& H, uint j,
				 uint numiter,
				 const mixsqp_control_params& control) {
  vec f = F.col(j);
  poismixsqp(L1,u,X.col(j),f,Z,H,numiter,control);
  return f;
}

// Perform one or more alt-SQP updates for a single column of the k x m
// factors matrix, F, when X is a sparse matrix.
inline vec altsqp_update_factor_sparse (const sp_mat& X, const mat& F,
					const mat& L1, const vec& u, mat& H,
					uint j, uint numiter,
					const mixsqp_control_params& control) {
  vec  x = nonzeros(X.col(j));
  uint n = x.n_elem;
  vec  f = F.col(j);
  uvec i(n);
  getcolnonzeros(X,i,j);
  poismixsqp(L1,u,x,i,f,H,numiter,control);
  return f;
}

// CLASS DEFINITIONS
// -----------------
// TO DO.

// FUNCTION DEFINITIONS
// --------------------
// Perform one or more alternating SQP ("alt-SQP") updates for the
// factors matrix, F, in which the matrix X is approximated by
// L*F. Input "numiter" specifies the number of EM updates to perform.
//
// Note that, unlike most other functions implemented in this package,
// the factors matrix F should be a k x m matrix, where k is the
// number factors, or "topics". This is done for ease of
// implementation, and for speed, because the computation is performed
// on the m columns of F, and the matrix is stored columnwise.
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
					     double numiter,
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
    Fnew.col(j) = altsqp_update_factor_sparse(X,F,L1,u,H,j,numiter,ctrl);
  return Fnew;
}
