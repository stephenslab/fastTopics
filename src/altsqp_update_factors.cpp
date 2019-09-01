#include "misc.h"
#include "altsqp.h"
#include <RcppParallel.h>

using namespace arma;

// INLINE FUNCTION DEFINITIONS
// ---------------------------
// Performs a single factor update (updates one column of the F matrix).
inline void altsqp_update_factor (const mat& X, const mat& F, const mat& L,
				  mat& Fnew, const vec& xscol, const vec& ls,
				  uint j, double e, uint numem, uint numsqp,
				  const mixsqp_control_params& control) {
  vec x = F.col(j);

  // Get the mixsqp inputs: set B to L[i,], and set w to X[i,j],
  // where i is the vector of indices such that X[i,j] > 0.
  uvec i = find(X.col(j) > 0);
  vec  w = nonzeros(X.col(j));
  uint n = i.n_elem;
  uint k = F.n_rows;
  mat  B(n,k);
  B = L.rows(i);
    
  // Run an SQP update.
  altsqp_update_em_sqp(B,w,ls,xscol(j),x,e,numem,numsqp,control);
      
  // Store the updated factors.
  Fnew.col(j) = x;
}

// Performs a single factor update (updates one column of the F
// matrix) when X is sparse.
inline void altsqp_update_factor_sparse (const sp_mat& X, const mat& F,
					 const mat& L, mat& Fnew,
					 const vec& xscol, const vec& ls,
					 uint j, double e, uint numem,
					 uint numsqp,
					 const mixsqp_control_params& control){
  vec x = F.col(j);

  // Get the mixsqp inputs: set B to L[i,], and set w to X[i,j],
  // where i is the vector of indices such that X[i,j] > 0.
  vec  w = nonzeros(X.col(j));
  uint n = w.n_elem;
  uint k = F.n_rows;
  uvec i(n);
  mat  B(n,k);
  getcolnonzeros(X,i,j);
  B = L.rows(i);
    
  // Run an SQP update.
  altsqp_update_em_sqp(B,w,ls,xscol(j),x,e,numem,numsqp,control);
      
  // Store the updated factors.
  Fnew.col(j) = x;
}

// CLASS DEFINITIONS
// -----------------
// This class is used to implement multithreaded computation of the
// factor updates in altsqp_update_factors_rcpp_parallel.
//
// [[Rcpp::depends(RcppParallel)]]
struct FactorUpdater : public RcppParallel::Worker {
  const mat& X;
  const mat& F;
  const mat& L;
  const vec& xscol;
  const vec& ls;
  mat&   Fnew;
  uint   numem;
  uint   numsqp;
  double e;
  mixsqp_control_params control;

  // This is used to create a FactorUpdater object.
  FactorUpdater (const arma::mat& X, const arma::mat& F, const arma::mat& L,
		 const arma::vec& xscol, const arma::vec& ls,
		 arma::mat& Fnew, double e, uint numem, uint numsqp,
		 mixsqp_control_params control) :
    X(X), F(F), L(L), xscol(xscol), ls(ls), Fnew(Fnew), numem(numem),
    numsqp(numsqp), e(e), control(control) { };

  // This function updates the factors for a given range of columns.
  void operator() (std::size_t begin, std::size_t end) {
    for (uint j = begin; j < end; j++)
      altsqp_update_factor(X,F,L,Fnew,xscol,ls,j,e,numem,numsqp,control);
  }
};

// This is a faster implementation of the R function altsqp.update.factors.
//
// The inputs and outputs differ slightly from the R function: X, an n
// x m matrix, must not be sparse; F is the transpose of the argument
// for the R function---an k x m matrix, in which F[j,] is the jth
// factor; L remains the same (n x k matrix); xscol is also the same,
// and must be equal to colSums(X); ls must be equal to colSums(L); e
// is control$e; and "control" is the same as in the R function.
//
// The return value is the k x m matrix of updated factors.
// 
// [[Rcpp::export]]
arma::mat altsqp_update_factors_rcpp (const arma::mat& X,
				      const arma::mat& F,
				      const arma::mat& L,
				      const arma::vec& xscol,
				      const arma::vec& ls, double e,
				      double numem, double numsqp,
				      Rcpp::List control) {
  mixsqp_control_params ctrl = get_mixsqp_control_params(control);

  // Initialize the return value.
  uint k = F.n_rows;
  uint m = F.n_cols;
  mat Fnew(k,m);

  // Repeat for each column of X (equivalently, for each column of F).
  for (uint j = 0; j < m; j++)
    altsqp_update_factor(X,F,L,Fnew,xscol,ls,j,e,numem,numsqp,ctrl);

  return Fnew;
}

// This is the same as altsqp_update_factors_rcpp, except that Intel
// Threading Building Blocks (TBB) are used to update the factors in
// parallel.
//
// [[Rcpp::export]]
arma::mat altsqp_update_factors_rcpp_parallel (const arma::mat& X,
					       const arma::mat& F,
					       const arma::mat& L,
					       const arma::vec& xscol,
					       const arma::vec& ls,
					       double e, double numem,
					       double numsqp,
					       Rcpp::List control) {
  mixsqp_control_params ctrl = get_mixsqp_control_params(control);

  // Initialize the return value.
  uint k = F.n_rows;
  uint m = F.n_cols;
  mat Fnew(k,m);
  
  // Create the worker.
  FactorUpdater worker(X,F,L,xscol,ls,Fnew,e,(uint) numem,(uint) numsqp,ctrl);
     
  // Update the factors with parallelFor.
  parallelFor(0,m,worker);
  
  return Fnew;
}

// This is the same as altsqp_update_factors_rcpp, except that input
// argument X is a sparse matrix.
//
// [[Rcpp::export]]
arma::mat altsqp_update_factors_sparse_rcpp (const arma::sp_mat& X,
					     const arma::mat& F,
					     const arma::mat& L,
					     const arma::vec& xscol,
					     const arma::vec& ls,
					     double e, double numem,
					     double numsqp,
					     Rcpp::List control) {
  mixsqp_control_params ctrl = get_mixsqp_control_params(control);

  // Initialize the return value.
  uint k = F.n_rows;
  uint m = F.n_cols;
  mat Fnew(k,m);

  // Repeat for each column of X (equivalently, for each column of F).
  for (uint j = 0; j < m; j++)
    altsqp_update_factor_sparse(X,F,L,Fnew,xscol,ls,j,e,numem,numsqp,ctrl);

  return Fnew;
}
