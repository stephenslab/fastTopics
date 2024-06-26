#include <RcppParallel.h>
#include "misc.h"
#include "poismix.h"

using namespace arma;

// FUNCTION DECLARATIONS
// ---------------------
void scd_update_factors (const mat& A, const mat& W, mat& H, const vec& j,
			 unsigned int numiter, double e);

void scd_update_factors_sparse (const sp_mat& A, const mat& W, mat& H,
				const vec& j, unsigned int numiter, 
				double e);

// INLINE FUNCTION DEFINITIONS
// ---------------------------
// Perform one or more SCD updates for a single column of H, in which
// A is approximated by the matrix product W*H.
inline void scd_update_factor (const mat& A, const mat& W, mat& H,
			       unsigned int j, unsigned int numiter, 
			       double e) {
  H.col(j) = scd_kl_update(W,A.col(j),H.col(j),numiter,e);
}

// Perform one or more SCD updates for a single column of H, in which
// A is approximated by the matrix product W*H, and A is a sparse
// matrix. Here, "sumw" should contain the precomputed column sums of
// W; that is, sumw = colSums(W).
inline void scd_update_factor_sparse (const sp_mat& A, const mat& W,
				      const vec& sumw, mat& H, unsigned int j,
				      unsigned int numiter, double e) {
  vec          a = nonzeros(A.col(j));
  unsigned int n = a.n_elem;
  uvec         i(n);
  getcolnonzeros(A,i,j);
  H.col(j) = scd_kl_update(W.rows(i),sumw,a,H.col(j),numiter,e);
}

// CLASS DEFINITIONS
// -----------------
// This class is used to implement multithreaded computation of the
// loadings updates in scd_update_factors_parallel_rcpp.
struct scd_factor_updater : public RcppParallel::Worker {
  const mat&   A;
  const mat&   W;
  mat&         H;
  const vec&   j;
  unsigned int numiter;
  double       e;

  // This is used to create a scd_factor_updater object.
  scd_factor_updater (const mat& A, const mat& W, mat& H, const vec& j,
		      unsigned int numiter, double e) :
    A(A), W(W), H(H), j(j), numiter(numiter), e(e) { };

  // This function updates the factors for a given range of columns.
  void operator() (std::size_t begin, std::size_t end) {
    for (unsigned int i = begin; i < end; i++)
      scd_update_factor(A,W,H,j(i),numiter,e);
  }
};

// This class is used to implement multithreaded computation of the
// loadings updates in scd_update_factors_sparse_parallel_rcpp.
struct scd_factor_updater_sparse : public RcppParallel::Worker {
  const sp_mat& A;
  const mat&    W;
  vec           sumw;
  mat&          H;
  const vec&    j;
  unsigned int  numiter;
  double        e;

  // This is used to create a scd_factor_updater_sparse object.
  scd_factor_updater_sparse (const sp_mat& A, const mat& W, mat& H, 
			     const vec& j, unsigned int numiter, double e) :
    A(A), W(W), sumw(W.n_cols), H(H), j(j), numiter(numiter), e(e) {
    sumw = trans(sum(W,0));
  };
 
  // This function updates the factors for a given range of columns.
  void operator() (std::size_t begin, std::size_t end) {
    for (unsigned int i = begin; i < end; i++)
      scd_update_factor_sparse(A,W,sumw,H,j(i),numiter,e);
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
				   const arma::mat& H, const arma::vec& j,
				   unsigned int numiter, double e) {
  mat Hnew = H;
  scd_update_factors(A,W,Hnew,j,numiter,e);
  return Hnew;
}

// Perform sequential co-ordinate descent (SCD) updates for the
// factors matrix (H) when the data matrix (A) is sparse.
//
// [[Rcpp::export]]
arma::mat scd_update_factors_sparse_rcpp (const arma::sp_mat& A,
					  const arma::mat& W,
					  const arma::mat& H,
					  const arma::vec& j,
					  unsigned int numiter, 
					  double e) {
  mat Hnew = H;
  scd_update_factors_sparse(A,W,Hnew,j,numiter,e);
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
					    const arma::vec& j,
					    unsigned int numiter, 
					    double e) {
  mat Hnew = H;
  scd_factor_updater worker(A,W,Hnew,j,numiter,e);
  parallelFor(0,j.n_elem,worker);
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
						   const arma::vec& j,
						   unsigned int numiter, 
						   double e) {
  mat Hnew = H;
  scd_factor_updater_sparse worker(A,W,Hnew,j,numiter,e);
  parallelFor(0,j.n_elem,worker);
  return Hnew;
}

// Iterate the SCD updates over all columns of H, in which A is
// approximated by the matrix product W*H.
void scd_update_factors (const mat& A, const mat& W, mat& H, const vec& j,
			 unsigned int numiter, double e) {
  unsigned int n = j.n_elem;
  for (unsigned int i = 0; i < n; i++)
    scd_update_factor(A,W,H,j(i),numiter,e);
}

// This is the same as scd_update_factors, except that the count data are
// stored as a sparse matrix.
void scd_update_factors_sparse (const sp_mat& A, const mat& W, mat& H,
				const vec& j, unsigned int numiter, double e) {
  unsigned int n    = j.n_elem;
  vec          sumw = trans(sum(W,0));
  for (unsigned int i = 0; i < n; i++)
    scd_update_factor_sparse(A,W,sumw,H,j(i),numiter,e);
}
