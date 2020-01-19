#include <RcppParallel.h>
#include "misc.h"
#include "poismix.h"

using namespace Rcpp;
using namespace arma;

// FUNCTION DECLARATIONS
// ---------------------
void scd_update (const mat& A, const mat& W, mat& H, uint numiter, double e);

void scd_update_sparse (const sp_mat& A, const mat& W, mat& H,
			uint numiter, double e);

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
    for (uint j = begin; j < end; j++)
      H.col(j) = scd_kl_update(W,A.col(j),H.col(j),numiter,e);
  }
};

// This class is used to implement multithreaded computation of the
// loadings updates in scd_update_factors_sparse_parallel_rcpp.
struct scd_factor_updater_sparse : public RcppParallel::Worker {
  const sp_mat& A;
  const mat&    W;
  mat&          H;
  uint          numiter;
  double        e;

  // This is used to create a scd_factor_updater object.
  scd_factor_updater_sparse (const sp_mat& A, const mat& W, mat& H,
			     uint numiter, double e) :
    A(A), W(W), H(H), numiter(numiter), e(e) { };

  // This function updates the factors for a given range of columns.
  void operator() (std::size_t begin, std::size_t end) {
    for (uint j = begin; j < end; j++) {
      vec  a = nonzeros(A.col(j));
      uint n = a.n_elem;
      uvec i(n);
      getcolnonzeros(A,i,j);
      H.col(j) = scd_kl_update_sparse(W,a,i,H.col(j),numiter,e);
    }
  }
};

// FUNCTION DEFINITIONS
// --------------------
// Perform sequential co-ordinate descent (SCD) updates for the
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

// Perform sequential co-ordinate descent (SCD) updates for the
// factors matrix (H) when the data matrix (A) is sparse.
//
// [[Rcpp::export]]
arma::mat scd_update_factors_sparse_rcpp (const arma::sp_mat& A,
					  const arma::mat& W,
					  const arma::mat& H,
					  uint numiter, double e) {
  mat Hnew = H;
  scd_update_sparse(A,W,Hnew,numiter,e);
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

// This does the same thing as scd_update_factors_sparse_rcpp, except
// that Intel Threading Building Blocks (TBB) are used to update the
// loadings in parallel.
//
// [[Rcpp::export]]
arma::mat scd_update_factors_sparse_parallel_rcpp (const arma::sp_mat& A, 
						   const arma::mat& W,
						   const arma::mat& H, 
						   uint numiter, double e) {
  mat Hnew = H;
  scd_factor_updater_sparse worker(A,W,Hnew,numiter,e);
  parallelFor(0,H.n_cols,worker);
  return Hnew;
}

// Iterate the SCD updates over all columns of H, in which A is
// approximated by the matrix product W*H.
void scd_update (const mat& A, const mat& W, mat& H, uint numiter, double e) {
  uint m = H.n_cols;
  for (uint j = 0; j < m; j++)
    H.col(j) = scd_kl_update(W,A.col(j),H.col(j),numiter,e);
}

// This is the same as scd_update, except that the count data are
// stored as a sparse matrix.
void scd_update_sparse (const sp_mat& A, const mat& W, mat& H,
			uint numiter, double e) {
  uint m = H.n_cols;
  for (uint j = 0; j < m; j++) {
    vec  a = nonzeros(A.col(j));
    uint n = a.n_elem;
    uvec i(n);
    getcolnonzeros(A,i,j);
    H.col(j) = scd_kl_update_sparse(W,a,i,H.col(j),numiter,e);
  }
}

