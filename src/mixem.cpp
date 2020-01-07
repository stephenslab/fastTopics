#include "mixem.h"
#include "misc.h"

using namespace arma;

// FUNCTION DEFINITIONS
// --------------------
// This is mainly used for testing the mixem C++ function.
//
// [[Rcpp::depends(RcppArmadillo)]]
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
//
// Note that x and L need not be normalized; they will automatically
// be normalized inside this function.
//
// Also note that it does not make sense to compute a MLE of the
// mixture proportions when n < 2 and/or when m < 2; mixem will supply
// a result in such cases, but the result will not be valid.
vec mixem (const mat& L, const vec& w, const vec& x0, uint numiter) {
  mat L1 = L;
  mat P  = L;
  vec x  = x0/sum(x0);
  normalizecols(L1);
  mixem(L1,w,x,P,numiter);
  return x;
}

// Use this variant of mixem if you plan on using the same L matrix
// multiple times, or for calling mixem multiple times with matrices
// of the same dimension. In the first case, you can reuse the L and P
// matrices; in the latter case, you can reuse the P matrix.
//
// For the result to be valid, the matrix L should be normalized
// beforehand so that each column sums to 1; P should be a matrix of
// the same size as L; and vector x should be normalized beforehand
// so that the entries sum to 1.
void mixem (const mat& L, const vec& w, vec& x, mat& P, uint numiter) {
  for (uint i = 0; i < numiter; i++)
    mixem_update(L,w,x,P);
}

// Perform a single EM update. See the mixem function for
// an explanation of the inputs.
void mixem_update (const mat& L, const vec& w, vec& x, mat& P) {
  double e = 1e-15;

  // Compute the posterior mixture assignment probabilities. A small
  // number is added to the posterior probabilities to prevent any
  // divisions by zero. This is the "E step".
  P = L;
  scalecols(P,x);
  normalizerowsbymax(P);
  P += e;
  normalizerows(P);
    
  // Update the mixture weights. This is the "M step".
  x = (trans(P) * w) / sum(w);
}
