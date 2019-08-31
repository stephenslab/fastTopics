#include "misc.h"
#include "altsqp.h"
#include <RcppParallel.h>

using namespace arma;

// CLASS DEFINITIONS
// -----------------
// TO DO: Explain here what this struct is for.
//
// [[Rcpp::depends(RcppParallel)]]
struct FactorUpdater : public RcppParallel::Worker {
  const mat& X;
  const mat& F;
  const mat& L;
  const vec& xscol;
  const vec& ls;
  mat&   Fnew;
  double e;
  uint   numem;
  uint   numsqp;
  mixsqp_control_params control;

  // This is used to create a FactorUpdater object.
  FactorUpdater (const arma::mat& X, const arma::mat& F, const arma::mat& L,
		 const arma::vec& xscol, const arma::vec& ls,
		 arma::mat& Fnew, double e, uint numem, uint numsqp,
		 mixsqp_control_params control) :
    X(X), F(F), L(L), xscol(xscol), ls(ls), Fnew(Fnew), e(e), numem(numem),
    numsqp(numsqp), control(control) { };

  // This function updates the factors for a given range of columns.
  void operator() (std::size_t begin, std::size_t end) {
    uint k = F.n_rows;
    for (std::size_t j = begin; j < end; j++) {
      vec x = F.col(j);

      // Get the mixsqp inputs: set B to L[i,], and set w to X[i,j],
      // where i is the vector of indices such that X[i,j] > 0.
      uvec i = find(X.col(j) > 0);
      vec  w = nonzeros(X.col(j));
      uint n = i.n_elem;
      mat  B(n,k);
      B = L.rows(i);
    
      // Run an SQP update.
      altsqp_update_em_sqp(B,w,ls,xscol(j),x,e,numem,numsqp,control);
      
      // Store the updated factors.
      Fnew.col(j) = x;
    }
  }
};

// FUNCTION DEFINITIONS
// --------------------
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
  for (uint j = 0; j < m; j++) {
    vec x = F.col(j);

    // Get the mixsqp inputs: set B to L[i,], and set w to X[i,j],
    // where i is the vector of indices such that X[i,j] > 0.
    uvec i = find(X.col(j) > 0);
    vec  w = nonzeros(X.col(j));
    uint n = i.n_elem;
    mat  B(n,k);
    B = L.rows(i);
    
    // Run an SQP update.
    altsqp_update_em_sqp(B,w,ls,xscol(j),x,e,(uint) numem,(uint) numsqp,ctrl);
      
    // Store the updated factors.
    Fnew.col(j) = x;
  }

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
  for (uint j = 0; j < m; j++) {
    vec x = F.col(j);

    // Get the mixsqp inputs: set B to L[i,], and set w to X[i,j],
    // where i is the vector of indices such that X[i,j] > 0.
    vec  w = nonzeros(X.col(j));
    uint n = w.n_elem;
    uvec i(n);
    mat B(n,k);
    getcolnonzeros(X,i,j);
    B = L.rows(i);
    
    // Run an SQP update.
    altsqp_update_em_sqp(B,w,ls,xscol(j),x,e,(uint) numem,(uint) numsqp,ctrl);
      
    // Store the updated factors.
    Fnew.col(j) = x;
  }

  return Fnew;
}

