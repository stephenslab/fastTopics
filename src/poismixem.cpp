#include "poismixem.h"
#include "misc.h"
#include "mixem.h"

using namespace arma;

// FUNCTION DEFINITIONS
// --------------------
// This is mainly used for testing the poismixem C++ function.
//
// [[Rcpp::export]]
arma::vec poismixem_rcpp (const arma::mat& L, const arma::vec& w,
		      const arma::vec& x0, uint numiter) {
  return poismixem(L,w,x0,numiter);
}

// Recover the mixture proportions of the multinomial mixture model
// from the mixture weights of the Poisson mixture model. The updated
// mixture proportions are stored in vector x. Additionally, upon
// completion, vector u contains column sums u[i] = sum(L[,i]), and
// matrix L is normalized so that each column sums to 1.
void p2m (mat& L, vec& u, vec& x) {
  u = sum(L,0);
  normalizecols(L);
  x %= u;
  x /= sum(x);
}

// Recover the mixture weights of the Poisson mixture model from the
// mixture weights of the multinomial mixture model, plus the "scale
// factor" (s). Input u should contain the column sums u[i] =
// sum(L[,i]).
void m2p (double s, const vec& u, vec& x) {
  x *= s;
  x /= u;
}

// Compute a maximum-likelihood estimate (MLE) of the mixture weights
// in a Poisson mixture model by iterating the multinomial mixture
// model EM updates for a fixed number of iterations.
//
// Input argument L is an n x m matrix with non-negative entries;
// input w is a vector of length n containing a count or "pseudocount"
// associated with each row of L; input argument x0 is the initial
// estimate of the mixture weights; input P is a matrix of the same
// dimension as L, and is used to store the posterior mixture
// assignment probabilities computed in the E-step; and input
// "numiter" specifies the number of EM updates to perform.
//
// Note that the second variant of "poismixem" modifies the L matrix;
// in particular, it normalizes the columns of L so that each column
// sums to 1.
//
// The return value is a vector of length m containing the updated
// mixture weights.
vec poismixem (const mat& L, const vec& w, const vec& x0, uint numiter) {
  mat L1 = L;
  mat P  = L;
  return poismixem(L1,w,x0,P,numiter);
}

vec poismixem (mat& L, const vec& w, const vec& x0, mat& P, uint numiter) {
  double s = sum(w);
  vec    x = x0;
  vec    u = x0;
  p2m(L,u,x);
  x = mixem(L,w,x,P,numiter);
  m2p(s,u,x);
  return x;
}
