#include <RcppParallel.h>
#include "misc.h"
#include "altsqp.h"

using namespace arma;

// INLINE FUNCTION DEFINITIONS
// ---------------------------
// Performs a single loadings update (updates one column of the L matrix).
inline void altsqp_update_loading (const mat& X, const mat& F, const mat& L,
				   mat& Lnew, const vec& xsrow, const vec& fs,
				   uint i, double e, uint numem, uint numsqp,
				   const mixsqp_control_params& control) {
  vec x = L.col(i);

  // Get the mixsqp inputs: set B to F[j,], and set w to X[j,i],
  // where i is the vector of indices such that X[j,i] > 0.
  uvec j = find(X.col(i) > 0);
  vec  w = nonzeros(X.col(i));
  uint m = j.n_elem;
  uint k = L.n_rows;
  mat  B(m,k);
  B = F.rows(j);

  // Run the EM and SQP updates.
  altsqp_update_em_sqp(B,w,fs,xsrow(i),x,e,numem,numsqp,control);

  // Store the updated loadings.
  Lnew.col(i) = x;
}

// Performs a single loadings update (updates one column of the L
// matrix) when X is sparse.
inline void altsqp_update_loading_sparse (const sp_mat& X, const mat& F,
					  const mat& L, mat& Lnew,
					  const vec& xsrow, const vec& fs,
					  uint i, double e, uint numem,
					  uint numsqp,
					  const mixsqp_control_params&control){
  vec x = L.col(i);

  // Get the mixsqp inputs: set B to F[j,], and set w to X[j,i],
  // where j is the vector of indices such that X[j,i] > 0.
  vec  w = nonzeros(X.col(i));
  uint m = w.n_elem;
  uint k = L.n_rows;
  uvec j(m);
  mat B(m,k);
  getcolnonzeros(X,j,i);
  B = F.rows(j);

  // Run the EM and SQP updates.
  altsqp_update_em_sqp(B,w,fs,xsrow(i),x,e,numem,numsqp,control);

  // Store the updated loadings.
  Lnew.col(i) = x;
}

// CLASS DEFINITIONS
// -----------------
// This class is used to implement multithreaded computation of the
// loadings updates in altsqp_update_loadings_rcpp_parallel.
//
// [[Rcpp::depends(RcppParallel)]]
struct LoadingsUpdater : public RcppParallel::Worker {
  const mat& X;
  const mat& F;
  const mat& L;
  const vec& xsrow;
  const vec& fs;
  mat&   Lnew;
  uint   numem;
  uint   numsqp;
  double e;
  mixsqp_control_params control;

  // This is used to create a LoadingsUpdater object.
  LoadingsUpdater (const arma::mat& X, const arma::mat& F, const arma::mat& L,
		   const arma::vec& xsrow, const arma::vec& fs,
		   arma::mat& Lnew, double e, uint numem, uint numsqp,
		   mixsqp_control_params control) :
    X(X), F(F), L(L), xsrow(xsrow), fs(fs), Lnew(Lnew), numem(numem),
    numsqp(numsqp), e(e), control(control) { };

  // This function updates the loadings for a given range of columns.
  void operator() (std::size_t begin, std::size_t end) {
    for (uint i = begin; i < end; i++)
      altsqp_update_loading(X,F,L,Lnew,xsrow,fs,i,e,numem,numsqp,control);
  }
};

// This class is used to implement multithreaded computation of the
// loadings updates in altsqp_update_loadings_rcpp_parallel_sparse.
//
// [[Rcpp::depends(RcppParallel)]]
struct LoadingsUpdaterSparse : public RcppParallel::Worker {
  const sp_mat& X;
  const mat& F;
  const mat& L;
  const vec& xsrow;
  const vec& fs;
  mat&   Lnew;
  uint   numem;
  uint   numsqp;
  double e;
  mixsqp_control_params control;

  // This is used to create a LoadingsUpdateSparse object.
  LoadingsUpdaterSparse (const arma::sp_mat& X, const arma::mat& F,
			 const arma::mat& L, const arma::vec& xsrow,
			 const arma::vec& fs, arma::mat& Lnew, double e,
			 uint numem, uint numsqp,
			 mixsqp_control_params control) :
    X(X), F(F), L(L), xsrow(xsrow), fs(fs), Lnew(Lnew), numem(numem),
    numsqp(numsqp), e(e), control(control) { };

  // This function updates the loadings for a given range of columns.
  void operator() (std::size_t begin, std::size_t end) {
    for (uint i = begin; i < end; i++)
      altsqp_update_loading_sparse(X,F,L,Lnew,xsrow,fs,i,e,numem,numsqp,
				   control);
  }
};

// FUNCTION DEFINITIONS
// --------------------
// This is a faster implementation of the R function altsqp.update.loadings.
//
// The inputs and outputs differ slightly from the R function: X, an m
// x n matrix, must not be sparse; F remains the same (m x k matrix);
// L is the transpose of the argument for the R function---an k x n
// matrix, in which L[j,] is the set of loadings corresponding to the
// jth factor; xsrow is the same, and must be equal to rowSums(X); fs
// must be equal to colSums(F); e is control$e; and "control" is the
// same as in the R function.
//
// Importantly, the input X for this function is the *transpose* of
// the X inputted to altsqp_update_factors_rcpp.
// 
// The return value is the k x n matrix of updated loadings.
// 
// [[Rcpp::export]]
arma::mat altsqp_update_loadings_rcpp (const arma::mat& X,
				       const arma::mat& F,
				       const arma::mat& L,
				       const arma::vec& xsrow,
				       const arma::vec& fs, double e,
				       double numem, double numsqp,
				       Rcpp::List control) {
  mixsqp_control_params ctrl = get_mixsqp_control_params(control);

  // Initialize the return value.
  uint k = L.n_rows;
  uint n = L.n_cols;
  mat Lnew(k,n);

  // Repeat for each row of X (equivalently, for each column of L).
  for (uint i = 0; i < n; i++)
    altsqp_update_loading(X,F,L,Lnew,xsrow,fs,i,e,numem,numsqp,ctrl);

  return Lnew;
}

// This is the same as altsqp_update_loadings_rcpp, except that input
// argument X is a sparse matrix.
//
// [[Rcpp::export]]
arma::mat altsqp_update_loadings_sparse_rcpp (const arma::sp_mat& X,
					      const arma::mat& F,
					      const arma::mat& L,
					      const arma::vec& xsrow,
					      const arma::vec& fs, double e,
					      double numem, double numsqp,
					      Rcpp::List control) {
  mixsqp_control_params ctrl = get_mixsqp_control_params(control);

  // Initialize the return value.
  uint k = L.n_rows;
  uint n = L.n_cols;
  mat Lnew(k,n);

  // Repeat for each row of X (equivalently, for each column of L).
  for (uint i = 0; i < n; i++)
    altsqp_update_loading_sparse(X,F,L,Lnew,xsrow,fs,i,e,numem,numsqp,ctrl);

  return Lnew;
}

// This is the same as altsqp_update_loadings_rcpp, except that Intel
// Threading Building Blocks (TBB) are used to update the loadings in
// parallel.
//
// [[Rcpp::export]]
arma::mat altsqp_update_loadings_rcpp_parallel (const arma::mat& X,
						const arma::mat& F,
						const arma::mat& L,
						const arma::vec& xsrow,
						const arma::vec& fs,
						double e, double numem,
						double numsqp,
						Rcpp::List control) {
  mixsqp_control_params ctrl = get_mixsqp_control_params(control);

  // Initialize the return value.
  uint k = L.n_rows;
  uint n = L.n_cols;
  mat Lnew(k,n);
  
  // Create the worker.
  LoadingsUpdater worker(X,F,L,xsrow,fs,Lnew,e,(uint)numem,(uint)numsqp,ctrl);
     
  // Update the loadings with parallelFor.
  parallelFor(0,n,worker);
  
  return Lnew;
}

// This is the same as altsqp_update_loadings_rcpp_sparse, except that
// Intel Threading Building Blocks (TBB) are used to update the
// loadings in parallel.
//
// [[Rcpp::export]]
arma::mat altsqp_update_loadings_rcpp_parallel_sparse (const arma::sp_mat& X,
						       const arma::mat& F,
						       const arma::mat& L,
						       const arma::vec& xsrow,
						       const arma::vec& fs,
						       double e, double numem,
						       double numsqp,
						       Rcpp::List control) {
  mixsqp_control_params ctrl = get_mixsqp_control_params(control);

  // Initialize the return value.
  uint k = L.n_rows;
  uint n = L.n_cols;
  mat Lnew(k,n);
  
  // Create the worker.
  LoadingsUpdaterSparse worker(X,F,L,xsrow,fs,Lnew,e,(uint) numem,
			       (uint) numsqp,ctrl);
     
  // Update the loadings with parallelFor.
  parallelFor(0,n,worker);
  
  return Lnew;
}
