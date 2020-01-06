#include "nnlm.h"

// Using version 0.4-3 of the NNLM package.

//[[Rcpp::export]]
Rcpp::List c_nnmf (const arma::mat& A, const unsigned int k, arma::mat W,
		   arma::mat H, uint max_iter, uint inner_max_iter) {
  
  for(uint i = 0; i < max_iter; i++) {
    
    // update W
    update(W, H, A.t(), inner_max_iter);
    
    // update H
    update(H, W, A, inner_max_iter);
  }
  
  return Rcpp::List::create(Rcpp::Named("W") = W,Rcpp::Named("H") = H);
}
