#include "cost.h"
#include "misc.h"

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

// This function is the same as cost_rcpp, except that X must be a
// sparse matrix.
//
// [[Rcpp::export]]
double cost_sparse_rcpp (const arma::sp_mat& X, const arma::mat& A,
			 const arma::mat& B, double e) {
  return cost_sparse(X,A,B,e);
}

// Helper function for cost_rcpp.
double cost (const arma::mat& X, const arma::mat& A,
	     const arma::mat& B, double e) {
  int    n = X.n_rows;
  int    m = X.n_cols;
  int    k = A.n_cols;
  double f = 0;

  // Quantities used in the computations below.
  vec x(n);
  vec y(n);
  vec b(k);
  
  // Repeat for each column of X.
  for (int j = 0; j < m; j++) {

    // Get the jth column of X and the jth column of B.
    copycol(X,j,x);
    copycol(B,j,b);

    // This is equivalent to the following R code:
    //
    //   sum(y - X[,j]*log(y + e))
    //
    // where 
    // 
    //   y = A %*% B[,j]
    //
    y  = A * b;
    f += sum(y);
    f -= sum(x % log(y + e));
  }
  
  return f;
}

// Helper function for cost_sparse_rcpp.
double cost_sparse (const arma::sp_mat& X, const arma::mat& A,
		    const arma::mat& B, double e) {
  int    n = X.n_rows;
  int    m = X.n_cols;
  int    k = A.n_cols;
  double f = 0;
  
  // Quantities used in the computations below.
  vec x(n);
  vec y(n);
  vec b(k);
  
  // Repeat for each column of X.
  for (int j = 0; j < m; j++) {

    // Initialize an iterator for the nonzero elements in the jth
    // column of X.
    sp_mat::const_col_iterator x   = X.begin_col(j);
    sp_mat::const_col_iterator end = X.end_col(j);

    // Get the jth column of B.
    copycol(B,j,b);

    // This is equivalent to the following R code:
    //
    //   sum(y - X[,j]*log(y + e))
    //
    // where 
    // 
    //   y = A %*% B[,j]
    //
    y  = A * b;
    f += sum(y);
    for(; x != end; ++x)
      f -= (*x) * log(y(x.row()) + e);
  }
  
  return f;
}
