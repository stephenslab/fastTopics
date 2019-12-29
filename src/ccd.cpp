#include "misc.h"

using namespace Rcpp;

// FUNCTION DECLARATIONS
// ---------------------
void ccd_update_factor (uint m, uint k, double* wt, double* wht,
			const double* vt, const double* H, double e);

void ccd_update_factors (uint n, uint m, uint k, const double* V,
			 const double* W, double* H, double* WH, double e);

void ccd (uint n, uint m, uint k, const double* V, double* W, double* H, 
	  double* WH, double* vt, double* wht, double e);

// FUNCTION DEFINITIONS
// --------------------
// [[Rcpp::export]]
void ccd_rcpp (const NumericMatrix& V, NumericMatrix& W, NumericMatrix& H,
	       NumericMatrix& WH, NumericVector& wht, NumericVector& vt,
	       double e) {
  ccd(V.nrow(),V.ncol(),W.nrow(),V.begin(),W.begin(),H.begin(),WH.begin(),
      vt.begin(),wht.begin(),e);
}

// [[Rcpp::export]]
void ccd_update_factors_rcpp (const NumericMatrix& V, const NumericMatrix& W, 
			      NumericMatrix& H, NumericMatrix& WH, double e) {
  ccd_update_factors(V.nrow(),V.ncol(),W.nrow(),V.begin(),W.begin(),
		     H.begin(),WH.begin(),e);
}

// Implements the cyclic co-ordinate descent updates described in
// Hsieh & Dhillon (2011).
void ccd (uint n, uint m, uint k, const double* V, double* W, double* H, 
	  double* WH, double* vt, double* wht, double e) {

  // Update W.
  for (uint i = 0; i < n; i++) {
    for (uint j = 0; j < m; j++) {
      wht[j] = WH[j*n + i];
      vt[j]  = V[j*n + i];
    }
    ccd_update_factor(m,k,W + i*k,wht,vt,H,e);
    for (uint j = 0; j < m; j++)
      WH[j*n + i] = wht[j];
  }

  // Update H.
  for (uint i = 0; i < m; i++)
    ccd_update_factor(n,k,H + i*k,WH + i*n,V + i*n,W,e);
}

// Implements the cyclic co-ordinate descent (CCD) update for the
// loadings matrix (W).
void ccd_update_factors (uint n, uint m, uint k, const double* V,
			 const double* W, double* H, double* WH, double e) {
  for (uint j = 0; j < m; j++)
    ccd_update_factor(n,k,H + j*k,WH + j*n,V + j*n,W,e);
}

// Implements the core part of the cyclic co-ordinate descent (CCD)
// updates. This implementation is adapted from the C++ code developed
// by Cho-Jui Hsieh and Inderjit Dhillon, and available for download
// at www.cs.utexas.edu/~cjhsieh/nmf.
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
