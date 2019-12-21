#include "mixem.h"
#include "misc.h"

using namespace arma;

// FUNCTION DEFINITIONS
// --------------------

// This is mainly used to test the mixem C++ function.
//
// [[Rcpp::export]]
arma::vec mixem_rcpp (const arma::mat& L, const arma::vec& w,
		      const arma::vec& x0, uint numiter, double e) {
  return mixem(L,w,x0,numiter,e);
}

// Compute maximum-likelihood estimates of the mixture proportions in
// a mixture model by iterating the EM updates for a fixed number of
// iterations.
vec mixem (const mat& L, const vec& w, const vec& x0, uint numiter, double e) {
  vec  x = x0;
  uint n = L.n_rows;
  uint m = L.n_cols;
  mat  P(n,m);
  
  // Iterate the E and M steps.
  for (uint i = 0; i < numiter; i++)
    mixem_update(L,w,x,P,e);
  
  return x;
}

// Perform a single EM update.
void mixem_update (const mat& L, const vec& w, vec& x, mat& P, double e) {

  // Compute the posterior mixture assignment probabilities. This is
  // the "E step".
  P = L;
  scalecols(P,x);
  P = P + e;
  normalizerows(P);
  
  // Update the mixture weights. This is the "M step".
  x = (trans(P) * w) / sum(w);
}
