#include <RcppParallel.h>
#include "misc.h"
#include "poismix.h"

using namespace Rcpp;
using namespace arma;

// FUNCTION DECLARATIONS
// ---------------------
void ccd_update_factors (const mat& V, const mat& W, mat& H, double e);

void ccd_update_factors_sparse (const sp_mat& V, const mat& W, mat& H, 
				double e);

// INLINE FUNCTION DEFINITIONS
// ---------------------------
// Perform a CCD update for a single column of H, in which V is
// approximated by the matrix product W*H.
inline void ccd_update_factor (const mat& V, const mat& W, mat& H,
			       uint j, double e) {
  H.col(j) = ccd_kl_update(W,V.col(j),H.col(j),e);
}

// Perform a SCD update for a single column of H, in which V is
// approximated by the matrix product W*H, and V is a sparse matrix.
// Here, "sumw" should contain the precomputed column sums of W; that
// is, sumw = colSums(W).
inline void ccd_update_factor_sparse (const sp_mat& V, const mat& W,
				      const vec& sumw, mat& H, uint j,
				      double e) {
  vec  v = nonzeros(V.col(j));
  uint n = v.n_elem;
  uvec i(n);
  getcolnonzeros(V,i,j);
  H.col(j) = ccd_kl_update(W.rows(i),sumw,v,H.col(j),e);
}

// FUNCTION DEFINITIONS
// --------------------
// Perform a single cyclic co-ordinate descent (CCD) update for the
// factors matrix (H). 
//
// The inputs are: V, the n x m data matrix; W, the n x k loadings
// matrix; H, the initial estimate of the k x m factors matrix; and
// e, a non-negative scalar specifying the minimum value of the
// updated parameters.
//
// This implementation is adapted from the C++ code developed by
// Cho-Jui Hsieh and Inderjit Dhillon, which is available for download
// at www.cs.utexas.edu/~cjhsieh/nmf.
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat ccd_update_factors_rcpp (const arma::mat& V, const arma::mat& W, 
				   arma::mat& H, double e) {
  mat Hnew = H;
  ccd_update_factors(V,W,Hnew,e);
  return Hnew;
}

// Perform sequential co-ordinate descent (CCD) updates for the
// factors matrix (H) when the data matrix (W) is sparse.
//
// [[Rcpp::export]]
arma::mat ccd_update_factors_sparse_rcpp (const arma::sp_mat& V,
					  const arma::mat& W,
					  const arma::mat& H, double e) {
  mat Hnew = H;
  ccd_update_factors_sparse(V,W,Hnew,e);
  return Hnew;
}

// Iterate the CCD updates over all columns of H, in which V is
// approximated by the matrix product W*H.
void ccd_update_factors (const mat& V, const mat& W, mat& H, double e) {
  uint m = H.n_cols;
  for (uint j = 0; j < m; j++)
    ccd_update_factor(V,W,H,j,e);
}

// This is the same as ccd_update_factors, except that the count data are
// stored as a sparse matrix.
void ccd_update_factors_sparse (const sp_mat& V, const mat& W, mat& H,
				double e) {
  uint m    = H.n_cols;
  vec  sumw = sum(W,0);
  for (uint j = 0; j < m; j++)
    ccd_update_factor_sparse(V,W,sumw,H,j,e);
}
