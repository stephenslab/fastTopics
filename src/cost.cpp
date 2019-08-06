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

// This is the same as cost_rcpp, except that input matrix X must be
// sparse.
//
// [[Rcpp::export]]
double cost_sparse_rcpp (const arma::sp_mat& X, const arma::mat& A,
			 const arma::mat& B, double e) {
  return 0;
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
    copycol(X,j,x);
    copycol(B,j,b);
    y = A * b;
    f += sum(y - x % log(y + e));
  }
  
  return f;
}
