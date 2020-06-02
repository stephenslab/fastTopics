#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List binom_stats_rcpp (const arma::mat& X, const arma::mat& F,
		       const arma::mat& L) {
 
  // Get the number of rows (n) and columns (m) of the counts matrix,
  // and the number of topics (k).
  unsigned int n = X.n_rows;
  unsigned int m = F.n_rows;
  unsigned int k = F.n_cols;
  
  // Initialize the outputs.
  mat    a1(m,k,fill::zeros);
  mat    a0(m,k,fill::zeros);
  rowvec b1(k,fill::zeros);
  rowvec b0(k,fill::zeros);

  // These variables are used to store intermediate results.
  rowvec p(k);
  double x;
  
  // Repeat for each row and column of the counts matrix.
  for (unsigned int i = 0; i < n; i++) 
    for (unsigned int j = 0; j < m; j++) {
      x = X(i,j);
    
      // Compute the posterior topic assignment probabilities.
      p = L.row(i) % F.row(j);
      p /= sum(p);
        
      // Update the expectations.
      a0.row(j) += x* (1-p);
      a1.row(j) += x*p;
      b0        += x*(1-p);
      b1        += x*p;
    }

  // Output the expectations.
  return List::create(Named("a1") = a1,
                      Named("a0") = a0,
		      Named("b1") = b1,
		      Named("b0") = b0);
}
