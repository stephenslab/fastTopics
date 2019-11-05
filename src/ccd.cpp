#include <Rcpp.h>

using namespace Rcpp;

#define max(a,b) ((a) > (b) ? (a) : (b))

// FUNCTION DECLARATIONS
// ---------------------
double update (int m, int k, double* wt, double* wht, const double* vt, 
	       double* H, double e);
void   ccd    (int n, int m, int k, const double* V, double* W, double* H, 
	       double* WH, double* vt, double* wht, double e);

// FUNCTION DEFINITIONS
// --------------------
// [[Rcpp::export]]
void ccd_rcpp (const NumericMatrix& V, NumericMatrix& W, NumericMatrix& H,
	       NumericMatrix& WH, NumericVector& wht, NumericVector& vt,
	       double e) {

  // Get the dimensions of the data matrix (n x m), and the number of
  // factors in the non-negative matrix factorization (k).
  R_xlen_t n = V.nrow();
  R_xlen_t m = V.ncol();
  R_xlen_t k = W.nrow();

  // Run the CCD updates.
  ccd(n,m,k,V.begin(),W.begin(),H.begin(),WH.begin(),vt.begin(),
      wht.begin(),e);
}

// Implements the cyclic co-ordinate descent updates described in
// Hsieh & Dhillon (2011).
void ccd (int n, int m, int k, const double* V, double* W, double* H, 
	  double* WH, double* vt, double* wht, double e) {
  double d, dh, dw;

  // Update W.
  dw = 0;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      wht[j] = WH[j*n + i];
      vt[j]  = V[j*n + i];
    }
    d  = update(m,k,W + i*k,wht,vt,H,e);
    dw = max(dw,d);
    for (int j = 0; j < m; j++)
      WH[j*n + i] = wht[j];
  }

  // Update H.
  dh = 0;
  for (int i = 0; i < m; i++) {
    d  = update(n,k,H + i*k,WH + i*n,V + i*n,W,e);
    dh = max(dh,d);
  }
}

// Implements the co-ordinate descent updates for the factors (H) and
// loadings (W).
double update (int m, int k, double* wt, double* wht, const double* vt, 
	       double* H, double e) {
  int i, j, hi;
  double d, g, h, t, w0, w1;
  double dmax = 0;

  for (i = 0; i < k; i++) {
    g = 0;
    h = 0;
    for (j = 0, hi = i; j < m; j++, hi += k) {
      t  = vt[j]/(wht[j] + e);
      g += H[hi]*(1 - t);
      h += H[hi]*H[hi]*t/(wht[j] + e);
    }
    w0 = wt[i];
    w1 = wt[i] - g/h + e;
    if (w1 < e)
      w1 = e;
    d     = w1 - w0;
    dmax  = max(dmax,fabs(d));
    wt[i] = w1;
    for (j = 0; j < m; j++)
      wht[j] += d * H[j*k + i];
  }

  return dmax;
}
