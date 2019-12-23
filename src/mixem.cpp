#include "mixem.h"
#include "misc.h"

using namespace arma;

// FUNCTION DEFINITIONS
// --------------------
// This is mainly used for testing the mixem C++ function.
//
// [[Rcpp::export]]
arma::vec mixem_rcpp (const arma::mat& L, const arma::vec& w,
		      const arma::vec& x0, uint numiter) {
  return mixem(L,w,x0,numiter);
}

// Compute a maximum-likelihood estimate (MLE) of the mixture
// proportions in the multinomial mixture model by iterating the EM
// updates for a fixed number of iterations.
//
// Input argument L is an n x m matrix with non-negative entries;
// input w is a vector of length n containing a non-negative "weight"
// associated with each row of L; input argument x0 is the initial
// estimate of the mixture proportions; input P is a matrix of the
// same dimension as L, and is used to store the posterior mixture
// assignment probabilities; and input "numiter" specifies the number
// of EM updates to perform.
//
// The return value is a vector of length m containing the updated
// mixture proportions.
vec mixem (const mat& L, const vec& w, const vec& x0, uint numiter) {
  mat P = L;
  return mixem(L,w,x0,P,numiter);
}

vec mixem (const mat& L, const vec& w, const vec& x0, mat& P, uint numiter) {
  vec x = x0;
  for (uint i = 0; i < numiter; i++)
    mixem_update(L,w,x,P);
  return x;
}

// Perform a single EM update. See the comments attached to mixem for
// an explanation of the inputs.
void mixem_update (const mat& L, const vec& w, vec& x, mat& P) {
  double e = 1e-15;

  // Compute the posterior mixture assignment probabilities. This is
  // the "E step".
  P = L;
  scalecols(P,x);
  normalizerowsbymax(P);
  P += e;
  normalizerows(P);
    
  // Update the mixture weights. This is the "M step".
  x = (trans(P) * w) / sum(w);
}
