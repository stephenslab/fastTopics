#include <RcppArmadillo.h>
#include "cost.h"

using namespace arma;

// FUNCTION DEFINITIONS
// --------------------
// TO DO: Explain what this function does.
//
// [[Rcpp::export]]
double cost_rcpp (const arma::mat& X, const arma::mat& A,
		  const arma::mat& B, double e) {
  return cost(X,A,B,e);
}

// Helper function for cost_rcpp.
double cost (const arma::mat& X, const arma::mat& A,
	     const arma::mat& B, double e) {
  int    n = X.n_rows;
  int    m = X.n_cols;
  int    K = A.n_cols;
  double f = 0;
  double ab;
  
  // Repeat for each row and column of X.
  for (int j = 0; j < m; j++)
    for (int i = 0; i < n; i++) {
      ab = 0;
      for (int k = 0; k < K; k++)
	ab += A(i,k) * B(k,j);
      f += ab - X(i,j)*log(ab + e);
    }
  
  return f;
}
