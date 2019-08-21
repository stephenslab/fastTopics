#include "misc.h"
#include "cost.h"

using namespace arma;

// FUNCTION DEFINITIONS
// --------------------
// Compute the value of the cost function for non-negative matrix
// factorization, in which matrix X is approximated by the matrix
// product A*B. This is equivalent to the negative Poisson
// log-likelihood after removing terms that do not depend on A or B.
//
// [[Rcpp::export]]
double cost_rcpp (const arma::mat& X, const arma::mat& A,
		  const arma::mat& B, double e) {
  return cost(X,A,B,e);
}

// This is the same as cost_rcpp, except that X must be sparse.
//
// [[Rcpp::export]]
double cost_sparse_rcpp (const arma::sp_mat& X, const arma::mat& A,
			 const arma::mat& B, double e) {
  return cost_sparse(X,A,B,e);
}

// This is the helper function for cost_rcpp.
double cost (const mat& X, const mat& A, const mat& B, double e) {
  uint   n = X.n_rows;
  uint   m = X.n_cols;
  double f = 0;
  vec    y(n);
  
  // Repeat for each column of X.
  for (uint j = 0; j < m; j++) {

    // This is equivalent to the following code in R:
    //
    //   sum(y - X[,j]*log(y + e))
    //
    // where 
    // 
    //   y = A %*% B[,j]
    //
    y  = A * B.col(j);
    f += sum(y);
    f -= sum(X.col(j) % log(y + e));
  }
  
  return f;
}

// Helper function for cost_sparse_rcpp.
double cost_sparse (const sp_mat& X, const mat& A, const mat& B, double e) {
  int    n = X.n_rows;
  int    m = X.n_cols;
  double f = 0;
  vec    y(n);
  
  // Repeat for each column of X.
  for (int j = 0; j < m; j++) {

    // Initialize an iterator for the nonzero elements in the jth
    // column of X.
    sp_mat::const_col_iterator xj  = X.begin_col(j);
    sp_mat::const_col_iterator xm = X.end_col(j);

    // This is equivalent to the following code in R:
    //
    //   sum(y - X[,j]*log(y + e))
    //
    // where 
    // 
    //   y = A %*% B[,j]
    //
    y  = A * B.col(j);
    f += sum(y);
    for(; xj != xm; ++xj)
      f -= (*xj) * log(y(xj.row()) + e);
  }
  
  return f;
}
