#include "misc.h"

using namespace Rcpp;

// FUNCTION DECLARATIONS
// ---------------------
void ccd_update_factor (uint m, uint k, double* wt, double* wht,
			const double* vt, const double* H, double e);

void ccd_update_factors (uint n, uint m, uint k, const double* V,
			 const double* W, double* H, double* WH, double e);

// FUNCTION DEFINITIONS
// --------------------
// Perform a single cyclic co-ordinate descent (CCD) update for the
// factors matrix (H). 
//
// The inputs are: V, the n x m data matrix; W, the loadings matrix,
// stored as a k x n matrix; H, the initial estimate of the k x m
// factors matrix; WH = crossprod(W,H), an n x m matrix; and e, a
// non-negative scalar specifying the minimum value of the updated
// parameters. Note that inputs H and WH are both modified.
//
// This implementation is adapted from the C++ code developed by
// Cho-Jui Hsieh and Inderjit Dhillon, which is available for download
// at www.cs.utexas.edu/~cjhsieh/nmf.
//
// [[Rcpp::export]]
void ccd_update_factors_rcpp (const NumericMatrix& V, const NumericMatrix& W, 
			      NumericMatrix& H, NumericMatrix& WH, double e) {
  ccd_update_factors(V.nrow(),V.ncol(),W.nrow(),V.begin(),W.begin(),
		     H.begin(),WH.begin(),e);
}

// Iterate the CCD updates over all columns of H.
void ccd_update_factors (uint n, uint m, uint k, const double* V,
			 const double* W, double* H, double* WH, double e) {
  for (uint j = 0; j < m; j++)
    ccd_update_factor(n,k,H + j*k,WH + j*n,V + j*n,W,e);
}

// Implements the core part of the cyclic co-ordinate descent (CCD)
// updates.
void ccd_update_factor (uint m, uint k, double* wt, double* wht,
			const double* vt, const double* H, double e) {
  double d, g, h, t, w0, w1;
  for (uint i = 0; i < k; i++) {
    g = 0;
    h = 0;
    for (uint j = 0, hi = i; j < m; j++, hi += k) {
      t  = vt[j]/(wht[j] + e);
      g += H[hi]*(1 - t);
      h += H[hi]*H[hi]*t/(wht[j] + e);
    }
    w0    = wt[i];
    w1    = wt[i] - g/h + e;
    w1    = maximum(w1,e);
    d     = w1 - w0;
    wt[i] = w1;
    for (uint j = 0; j < m; j++)
      wht[j] += d * H[j*k + i];
  }
}
