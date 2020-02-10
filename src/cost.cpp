#include "misc.h"
#include "cost.h"

using namespace arma;

// FUNCTION DEFINITIONS
// --------------------
// Compute negative log-likelihoods for assessing a topic model fit or
// quality of a non-negative matrix factorization, in which matrix X
// is approximated by matrix product A * B.
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec cost_rcpp (const arma::mat& X, const arma::mat& A,
		     const arma::mat& B, double e, bool poisson) {
  return cost(X,A,B,e,poisson);
}

// This is the same as cost_rcpp, except that X must be sparse.
//
// [[Rcpp::export]]
arma::vec cost_sparse_rcpp (const arma::sp_mat& X, const arma::mat& A,
	  		    const arma::mat& B, double e, bool poisson) {
  return cost_sparse(X,A,B,e,poisson);
}

// This is the helper function for cost_rcpp.
arma::vec cost (const mat& X, const mat& A, const mat& B, double e, 
             bool poisson) {
  unsigned int n = X.n_rows;
  unsigned int m = X.n_cols;
  vec  f(n,fill::zeros);
  vec  y(n);
  
  // Repeat for each column of X.
  for (unsigned int j = 0; j < m; j++) {

    // This is equivalent to the following R code:
    //
    //   f = f + poisson*y - X[,j]*log(y + e))
    //
    // where 
    // 
    //   y = A %*% B[,j]
    //
    y  = A * B.col(j);
    f -= X.col(j) % log(y + e);
    if (poisson)
      f += y;
  }
  
  return f;
}

// Helper function for cost_sparse_rcpp.
arma::vec cost_sparse (const sp_mat& X, const mat& A, const mat& B,
		       double e, bool poisson) {
  unsigned int n = X.n_rows;
  unsigned int m = X.n_cols;
  unsigned int i;
  vec  f(n,fill::zeros);
  vec  y(n);
  
  // Repeat for each column of X.
  for (unsigned int j = 0; j < m; j++) {

    // Initialize an iterator for the nonzero elements in the jth
    // column of X.
    sp_mat::const_col_iterator xj = X.begin_col(j);
    sp_mat::const_col_iterator xm = X.end_col(j);

    // This is equivalent to the following R code:
    //
    //   f = f + poisson*y - X[,j]*log(y + e)
    //
    // where 
    // 
    //   y = A %*% B[,j]
    //
    y = A * B.col(j);
    for(; xj != xm; ++xj) {
      i     = xj.row();
      f(i) -= (*xj) * log(y(i) + e);
    }
    if (poisson)
      f += y;
  }
  
  return f;
}
