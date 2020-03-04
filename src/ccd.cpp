#include <RcppParallel.h>
#include "misc.h"
#include "poismix.h"

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
			       unsigned int j, double e) {
  H.col(j) = ccd_kl_update(W,V.col(j),H.col(j),e);
}

// Perform a CCD update for a single column of H, in which V is
// approximated by the matrix product W*H, and V is a sparse matrix.
// Here, "sumw" should contain the precomputed column sums of W; that
// is, sumw = colSums(W).
inline void ccd_update_factor_sparse (const sp_mat& V, const mat& W,
				      const vec& sumw, mat& H, 
				      unsigned int j, double e) {
  vec          v = nonzeros(V.col(j));
  unsigned int n = v.n_elem;
  uvec         i(n);
  getcolnonzeros(V,i,j);
  H.col(j) = ccd_kl_update(W.rows(i),sumw,v,H.col(j),e);
}

// CLASS DEFINITIONS
// -----------------
// This class is used to implement multithreaded computation of the
// loadings updates in ccd_update_factors_parallel_rcpp.
struct ccd_factor_updater : public RcppParallel::Worker {
  const mat& V;
  const mat& W;
  mat&       H;
  double     e;

  // This is used to create a ccd_factor_updater object.
  ccd_factor_updater (const mat& V, const mat& W, mat& H, double e) :
    V(V), W(W), H(H), e(e) { };

  // This function updates the factors for a given range of columns.
  void operator() (std::size_t begin, std::size_t end) {
    for (unsigned int j = begin; j < end; j++)
      ccd_update_factor(V,W,H,j,e);
  }
};

// This class is used to implement multithreaded computation of the
// loadings updates in ccd_update_factors_sparse_parallel_rcpp.
struct ccd_factor_updater_sparse : public RcppParallel::Worker {
  const sp_mat& V;
  const mat&    W;
  vec           sumw;
  mat&          H;
  double        e;

  // This is used to create a ccd_factor_updater_sparse object.
  ccd_factor_updater_sparse (const sp_mat& V, const mat& W, mat& H, double e) :
    V(V), W(W), sumw(W.n_cols), H(H), e(e) {
    sumw = sum(W,0);
  };

  // This function updates the factors for a given range of columns.
  void operator() (std::size_t begin, std::size_t end) {
    for (unsigned int j = begin; j < end; j++)
      ccd_update_factor_sparse(V,W,sumw,H,j,e);
  }
};

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

// This does the same thing as ccd_update_factors_rcpp, except that
// Intel Threading Building Blocks (TBB) are used to update the
// loadings in parallel.
//
// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::export]]
arma::mat ccd_update_factors_parallel_rcpp (const arma::mat& V, 
					    const arma::mat& W,
					    const arma::mat& H, double e) {
  mat Hnew = H;
  ccd_factor_updater worker(V,W,Hnew,e);
  parallelFor(0,H.n_cols,worker);
  return Hnew;
}

// This does the same thing as ccd_update_factors_sparse_rcpp, except
// that Intel Threading Building Blocks (TBB) are used to update the
// loadings in parallel.
//
// [[Rcpp::export]]
arma::mat ccd_update_factors_sparse_parallel_rcpp (const arma::sp_mat& V, 
						   const arma::mat& W,
						   const arma::mat& H, 
						   double e) {
  mat Hnew = H;
  ccd_factor_updater_sparse worker(V,W,Hnew,e);
  parallelFor(0,H.n_cols,worker);
  return Hnew;
}

// Iterate the CCD updates over all columns of H, in which V is
// approximated by the matrix product W*H.
void ccd_update_factors (const mat& V, const mat& W, mat& H, double e) {
  unsigned int m = H.n_cols;
  for (unsigned int j = 0; j < m; j++)
    ccd_update_factor(V,W,H,j,e);
}

// This is the same as ccd_update_factors, except that the count data are
// stored as a sparse matrix.
void ccd_update_factors_sparse (const sp_mat& V, const mat& W, mat& H,
				double e) {
  unsigned int m = H.n_cols;
  vec sumw = sum(W,0);
  for (unsigned int j = 0; j < m; j++)
    ccd_update_factor_sparse(V,W,sumw,H,j,e);
}
