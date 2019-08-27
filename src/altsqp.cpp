#include "misc.h"
#include "mixem.h"
#include "mixsqp.h"
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// FUNCTION DECLARATIONS
// ---------------------
vec  get_modified_problem_params (mat& B, vec& w, const vec& bs, double ws,
				  vec& x);
void altsqp_update_em_sqp (mat& B, vec& w, const vec& bs, double ws, vec& x,
			   double e, uint numem, uint numsqp,
			   mixsqp_control_params control);

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
				      const arma::vec& ls,
				      double e, double numem,
				      double numsqp, List control) {
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

    // This is also a good point to check for a user interrupt; if the
    // user requests an interrupt, then an exception is thrown and
    // control is returned to the R console.
    Rcpp::checkUserInterrupt();
  }

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
					     double numsqp, List control) {
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

    // This is also a good point to check for a user interrupt; if the
    // user requests an interrupt, then an exception is thrown and
    // control is returned to the R console.
    Rcpp::checkUserInterrupt();
  }

  return Fnew;
}

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
				       const arma::vec& fs,
				       double e, double numem,
				       double numsqp, List control) {
  mixsqp_control_params ctrl = get_mixsqp_control_params(control);

  // Initialize the return value.
  uint k = L.n_rows;
  uint n = L.n_cols;
  mat Lnew(k,n);

  // Repeat for each row of X (equivalently, for each column of L).
  for (uint i = 0; i < n; i++) {
    vec x = L.col(i);

    // Get the mixsqp inputs: set B to F[j,], and set w to X[j,i],
    // where i is the vector of indices such that X[j,i] > 0.
    uvec j = find(X.col(i) > 0);
    vec  w = nonzeros(X.col(i));
    uint m = j.n_elem;
    mat  B(m,k);
    B = F.rows(j);

    // Run one EM update and one SQP update.
    altsqp_update_em_sqp(B,w,fs,xsrow(i),x,e,(uint) numem,(uint) numsqp,ctrl);

    // Store the updated loadings.
    Lnew.col(i) = x;
    
    // This is also a good point to check for a user interrupt; if the
    // user requests an interrupt, then an exception is thrown and
    // control is returned to the R console.
    Rcpp::checkUserInterrupt();
  }

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
					      const arma::vec& fs,
					      double e, double numem,
					      double numsqp, List control) {
  mixsqp_control_params ctrl = get_mixsqp_control_params(control);

  // Initialize the return value.
  uint k = L.n_rows;
  uint n = L.n_cols;
  mat Lnew(k,n);

  // Repeat for each row of X (equivalently, for each column of L).
  for (uint i = 0; i < n; i++) {
    vec x = L.col(i);

    // Get the mixsqp inputs: set B to F[j,], and set w to X[j,i],
    // where j is the vector of indices such that X[j,i] > 0.
    vec  w = nonzeros(X.col(i));
    uint m = w.n_elem;
    uvec j(m);
    mat B(m,k);
    getcolnonzeros(X,j,i);
    B = F.rows(j);

    // Run one EM update and one SQP update.
    altsqp_update_em_sqp(B,w,fs,xsrow(i),x,e,(uint) numem,(uint) numsqp,ctrl);

    // Store the updated loadings.
    Lnew.col(i) = x;
    
    // This is also a good point to check for a user interrupt; if the
    // user requests an interrupt, then an exception is thrown and
    // control is returned to the R console.
    Rcpp::checkUserInterrupt();
  }

  return Lnew;
}

// Run one EM update and SQP update on the modified problem, then
// recover the updated solution to the unmodified problem.
void altsqp_update_em_sqp (mat& B, vec& w, const vec& bs, double ws, vec& x,
			   double e, uint numem, uint numsqp,
			   mixsqp_control_params control) {

  // Update the solution to the modified problem.
  uint n = B.n_rows;
  vec  y = get_modified_problem_params(B,w,bs,ws,x);
  vec  ev(n);
  ev.fill(e);
  if (numem > 0)
    mixem(B,w,x,ev,numem);
  if (numsqp > 0)
    mixsqp(B,w,x,ev,numsqp,control,false);

  // Recover the updated solution to the unmodified problem.
  x %= y;
}

// Set up the modified problem solved using either EM or the mix-SQP
// algorithm.
vec get_modified_problem_params (mat& B, vec& w, const vec& bs,
				 double ws, vec& x) {
  uint m = B.n_cols;
  vec y(m);
  y.fill(ws);
  y /= bs;
  w /= ws;
  x /= y;
  scalecols(B,y);
  return y;
}
