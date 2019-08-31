#include "misc.h"
#include "altsqp.h"
#include <RcppParallel.h>

using namespace arma;

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
					      const arma::vec& fs, double e,
					      double numem, double numsqp,
					      Rcpp::List control) {
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
  }

  return Lnew;
}
