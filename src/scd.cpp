#include "misc.h"

using namespace Rcpp;
using namespace arma;

// FUNCTION DECLARATIONS
// ---------------------
void scd_update (mat& H, const mat& Wt, const mat& A, uint numiter, double e);

void scd_kl_update (subview_col<double> Hj, const mat& Wt, const vec& Aj,
		    const vec& sumW, uint numiter, double e);

// FUNCTION DEFINITIONS
// --------------------
// Perform a sequential co-ordinate descent (SCD) update for the
// loadings matrix (W).
//
// The inputs are: A, the m x n matrix to be decomposed, in which A is
// approximated as crossprod(H,W); W, the initial estimate of the k x
// n loadings matrix; H, the k x m factors matrix; numiter, the number
// of inner-loop iterations to perform; and e, a non-negative scalar
// included in the computations to prevent NaNs due to division by zero.
//
// This implementation is adapted from the R and C++ code developed by
// Xihui Lin and Paul Boutros, which is available for download at
// https://github.com/linxihui/NNLM.
//
// [[Rcpp::export]]
arma::mat scd_update_loadings_rcpp (const arma::mat& A, const arma::mat& W,
				    const arma::mat& H, uint numiter,
				    double e) {
  mat Wnew = W;
  scd_update(Wnew,H,A,numiter,e);
  return Wnew;
}

// Iterate the SCD updates over all columns of H.
void scd_update (mat& H, const mat& Wt, const mat& A, uint numiter, double e) {
  uint m    = A.n_cols;
  vec  sumw = sum(Wt,1);
  for (uint j = 0; j < m; j++) 
    scd_kl_update(H.col(j),Wt,A.col(j),sumw,numiter,e);
}

// Implements the core part of the sequential co-ordinate descent
// (SCD) updates.
void scd_kl_update (subview_col<double> Hj, const mat& Wt, const vec& Aj,
		    const vec& sumw, uint numiter, double e) {
  uint   k   = Wt.n_rows;
  vec    Ajt = Wt.t()*Hj;
  vec    mu;
  double a;
  double b;
  double x;

  for (uint t = 0; t < numiter; t++)
    for (uint i = 0; i < k; i++) {
      mu = Wt.row(i).t()/(Ajt + e);
      a  = dot(Aj,square(mu));
      b  = dot(Aj,mu) - sumw(i) + a*Hj(i);
      x  = b/(a + e);
      x  = maximum(x,0);
      if (x != Hj(i)) {
	Ajt   += (x - Hj(i)) * Wt.row(i).t();
	Hj(i)  = x;
      }
    }
}
