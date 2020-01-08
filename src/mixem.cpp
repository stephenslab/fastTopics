#include "mixem.h"
#include "misc.h"

using namespace arma;

// FUNCTION DECLARATIONS
// ---------------------
void mixem_update (const arma::mat& L1, const arma::vec& w,
		   arma::vec& x, arma::mat& P);

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
// Note that x0 and L need not be normalized; they will automatically
// be normalized inside this function.
//
// Also note that it does not make sense to compute a MLE of the
// mixture proportions when n < 2 and/or when m < 2; mixem will supply
// a result in such cases, but the result will not be valid.
vec mixem (const mat& L, const vec& w, const vec& x0, uint numiter) {
  mat L1 = L;
  mat P  = L;
  vec x  = x0;
  normalizecols(L1);
  mixem(L1,w,x,P,numiter);
  return x;
}

// Use this variant of mixem if you plan on using the same L matrix
// multiple times, or for calling mixem multiple times with matrices
// of the same dimension. In the first case, you can reuse the L1 and
// P matrices; in the latter case, you can reuse the P matrix.
//
// For the result to be valid, the matrix L1 should be normalized
// beforehand so that each column sums to 1. P should be a matrix of
// the same size as L1.
//
// Note that x need not be normalized; it will automatically be
// normalized inside this function.
void mixem (const mat& L1, const vec& w, vec& x, mat& P, uint numiter) {
  for (uint i = 0; i < numiter; i++)
    mixem_update(L1,w,x,P);
}

// Perform a single EM update. For this update to be valid, the matrix
// L1 should be normalized beforehand so that each column sums to 1.
// Note that x need not be normalized; it will automatically be
// normalized inside this function.
void mixem_update (const mat& L1, const vec& w, vec& x, mat& P) {
  double e = 1e-15;

  // Normalize the "weights".
  vec w1 = w/sum(w);

  // Normalize the mixture proportions.
  x /= sum(x);
  
  // Compute the posterior mixture assignment probabilities. A small
  // number is added to the posterior probabilities to prevent any
  // divisions by zero. This is the "E step".
  P = L1;
  scalecols(P,x);
  normalizerowsbymax(P);
  P += e;
  normalizerows(P);
    
  // Update the mixture weights. This is the "M step".
  x = trans(P) * w1;
}
