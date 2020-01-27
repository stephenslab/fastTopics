#include <RcppParallel.h>
#include "misc.h"
#include "poismix.h"

using namespace Rcpp;
using namespace arma;

// FUNCTION DECLARATIONS
// ---------------------
void ccd_update_factors (const mat& V, const mat& W, mat& H, double e);

// INLINE FUNCTION DEFINITIONS
// ---------------------------
// Perform one or more SCD updates for a single column of H, in which
// A is approximated by the matrix product W*H.
inline void ccd_update_factor (const mat& V, const mat& W, mat& H,
			       uint j, double e) {
  H.col(j) = ccd_kl_update(W,V.col(j),H.col(j),e);
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

// Iterate the CCD updates over all columns of H, in which V is
// approximated by the matrix product W*H.
void ccd_update_factors (const mat& V, const mat& W, mat& H, double e) {
  uint m = H.n_cols;
  for (uint j = 0; j < m; j++)
    ccd_update_factor(V,W,H,j,e);
}
