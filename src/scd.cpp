#include <RcppParallel.h>
#include "misc.h"

using namespace Rcpp;
using namespace arma;

// FUNCTION DECLARATIONS
// ---------------------
void scd_update (const mat& A, const mat& W, mat& H, uint numiter, double e);

void scd_update_sparse (mat& H, const mat& Wt, const sp_mat& A, uint numiter,
			double e);

void scd_kl_update (const vec& a, const mat& W, vec& h,
		    uint numiter, double e);

// INLINE FUNCTION DEFINITIONS
// ---------------------------
// Compute the first-order (g) and second-order (s) partial
// derivatives with respect to factor element H[i,j], where a = A[,j]
// is the jth column of the data matrix A, w = W[,i] is the ith column
// of the loadings matrix W, h = H[,j] is the jth column of the
// factors matrix H, and wh = W*h. Input argument u is a vector
// storing an intermediate result of the same size as w.
inline void compute_grad (const vec& a, const vec& w, const vec& h,
			  const vec& wh, uint i, double& g, double& s,
			  vec& u, double e) {
  u = w / (wh + e);
  s = dot(a,square(u));
  g = dot(a,u) - sum(w) + s*h(i);
}

// Given the first-order (g) and second-order (s) partial derivatives,
// compute the updated value of H[i,j] inside the factors matrix H.
inline double project_iterate (double g, double s, double e) {
  double h = g/(s + e);
  h = maximum(h,0);
  return h;
}

// CLASS DEFINITIONS
// -----------------
// This class is used to implement multithreaded computation of the
// loadings updates in scd_update_factors_parallel_rcpp.
struct scd_factor_updater : public RcppParallel::Worker {
  const mat& A;
  const mat& W;
  mat&       H;
  uint       numiter;
  double     e;

  // This is used to create a scd_factor_updater object.
  scd_factor_updater (const mat& A, const mat& W, mat& H, uint numiter,
		      double e) :
    A(A), W(W), H(H), numiter(numiter), e(e) { };

  // This function updates the factors for a given range of columns.
  void operator() (std::size_t begin, std::size_t end) {
    uint k = H.n_rows;
    vec  h(k);
    for (uint j = begin; j < end; j++) {
      h = H.col(j);
      scd_kl_update(A.col(j),W,h,numiter,e);
      H.col(j) = h;
    }
  }
};

// FUNCTION DEFINITIONS
// --------------------
// Perform a sequential co-ordinate descent (SCD) update for the
// factors matrix (H).
//
// The inputs are: A, the m x n matrix to be decomposed, in which A is
// approximated as W*H; W, the initial estimate of the n x k loadings
// matrix; H, the k x m factors matrix; numiter, the number of
// inner-loop iterations to perform; and e, a non-negative scalar
// included in the computations to prevent NaNs due to division by
// zero.
//
// This implementation is adapted from the R and C++ code developed by
// Xihui Lin and Paul Boutros, which is available for download at
// https://github.com/linxihui/NNLM.
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat scd_update_factors_rcpp (const arma::mat& A, const arma::mat& W,
				   const arma::mat& H, uint numiter,
				   double e) {
  mat Hnew = H;
  scd_update(A,W,Hnew,numiter,e);
  return Hnew;
}

// TO DO: Explain here what this function does, and how to use it.
//
// [[Rcpp::export]]
arma::mat scd_update_factors_sparse_rcpp (const arma::sp_mat& A,
					  const arma::mat& W,
					  const arma::mat& H,
					  uint numiter, double e) {
  mat Hnew = H;
  scd_update_sparse(Hnew,W,A,numiter,e);
  return Hnew;
}

// This does the same thing as scd_update_factors_rcpp, except that
// Intel Threading Building Blocks (TBB) are used to update the
// loadings in parallel.
//
// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::export]]
arma::mat scd_update_factors_parallel_rcpp (const arma::mat& A, 
					    const arma::mat& W,
					    const arma::mat& H, 
					    uint numiter, double e) {
  mat Hnew = H;
  scd_factor_updater worker(A,W,Hnew,numiter,e);
  parallelFor(0,H.n_cols,worker);
  return Hnew;
}

// Iterate the SCD updates over all columns of H (equivalently,
// columns of A), in which A is approximated by the matrix product W*H.
void scd_update (const mat& A, const mat& W, mat& H, uint numiter, double e) {
  uint k = H.n_rows;
  uint m = H.n_cols;
  vec  h(k);
  for (uint j = 0; j < m; j++) {
    h = H.col(j);
    scd_kl_update(A.col(j),W,h,numiter,e);
    H.col(j) = h;
  }
}

// TO DO: Explain here what this function does, and how to use it.
void scd_update_sparse (mat& H, const mat& W, const sp_mat& A,
			uint numiter, double e) {
  // TO DO.
}

// NOTES:
//
//   + A is approximated by matrix product W*H.
//   + a = column of A.
//   + h = column of H.
//  
// Implements the core part of the sequential co-ordinate descent
// (SCD) updates, in which vector a (column of matrix A) is
// approximated by the matrix-vector product W*h, where h is a column
// of the factors matrix H.
void scd_kl_update (const vec& a, const mat& W, vec& h,
		    uint numiter, double e) {
  uint   n  = W.n_rows;
  uint   k  = W.n_cols;
  vec    wh = W * h;
  vec    w(n);
  vec    u(n);
  double s, g, h0, h1;
  for (uint iter = 0; iter < numiter; iter++)
    for (uint i = 0; i < k; i++) {
      w = W.col(i);
      compute_grad(a,w,h,wh,i,g,s,u,e);
      h0 = h(i);
      h1 = project_iterate(g,s,e);
      if (h1 != h0)
	wh += (h1 - h0) * w;
      h(i) = h1;
    }
}

// TO DO: Explain here what this function does, and how to use it.
void scd_kl_update_sparse (uint numiter, double e) {
  // TO DO.
}
