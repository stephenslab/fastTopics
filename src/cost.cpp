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
double cost_rcpp (const arma::mat& X, const arma::mat& A, const arma::mat& B, 
		  double e) {
  return cost(X,A,B,e);
}

// This is the helper function for cost_rcpp.
double cost (const mat& X, const mat& A, const mat& B, double e) {
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
