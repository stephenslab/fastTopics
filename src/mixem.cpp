#include "mixem.h"
#include "misc.h"

using namespace arma;

// FUNCTION DEFINITIONS
// --------------------

// This is mainly used to test the mixem C++ function.
//
// [[Rcpp::export]]
arma::vec mixem_rcpp (const arma::mat& L, const arma::vec& w,
		      const arma::vec& x0, uint numiter) {
  mat P = L;
  return mixem(L,w,x0,P,numiter);
}

// Compute a maximum-likelihood estimate (MLE) of the mixture
// proportions in the multinomial mixture model by iterating the EM
// updates for a fixed number of iterations.
// 
// Input argument P is a matrix of the same dimension as L, and is
// used to store the posterior mixture assignment probabilities.
//
vec mixem (const mat& L, const vec& w, const vec& x0, mat& P, uint numiter) {
  vec x = x0;
  for (uint i = 0; i < numiter; i++)
    mixem_update(L,w,x,P);
  return x;
}

// Perform a single EM update.
void mixem_update (const mat& L, const vec& w, vec& x, mat& P) {
  double e = 1e-15;

  // Compute the posterior mixture assignment probabilities. This is
  // the "E step".
  P = L;
  scalecols(P,x);
  P += e;
  normalizerows(P);
  
  // Update the mixture weights. This is the "M step".
  x = (trans(P) * w) / sum(w);
}
