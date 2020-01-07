#include "poismix.h"
#include "misc.h"
#include "mixem.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;

// FUNCTION DEFINITIONS
// --------------------
// This is mainly used to test the first variant of the poismixem C++
// function.
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec poismixem_rcpp (const arma::mat& L, const arma::vec& w,
			  const arma::vec& x0, uint numiter) {
  return poismixem(L,w,x0,numiter);
}

// This is mainly used to test the second variant of the poismixem C++
// function.
//
// [[Rcpp::export]]
arma::vec poismixem2_rcpp (const arma::mat& L1, const arma::vec& w,
			   const arma::vec& u, const arma::vec& x0,
			   uint numiter) {
  vec x = x0;
  mat P = L1;
  poismixem(L1,u,w,x,P,numiter);
  return x;
}

// This is mainly used to test the third variant of the poismixem C++
// function.
//
// [[Rcpp::export]]
arma::vec poismixem3_rcpp (const arma::mat& L1, const arma::vec& w,
			   const arma::vec& u, const arma::uvec& i,
			   const arma::vec& x0, uint numiter) {
  vec x = x0;
  poismixem(L1,u,w,i,x,numiter);
  return x;
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
  vec u  = sum(L,0);
  vec x  = x0;
  normalizecols(L1);
  poismixem(L1,u,w,x,P,numiter);
  return x;
}

// Use this variant of poismixem if you plan on calling poismixem
// multiple times with the same matrix L. In this call, input u should
// contain the column sums, u = colSums(L), in which L is the matrix
// prior to normalization, and matrix L1 is the normalized version of
// L in which each column sums to 1; that is, L1 = normalize.cols(L).
void poismixem (const mat& L1, const vec& u, const vec& w, vec& x, mat& P, 
		uint numiter) {

  // Recover the mixture proportions of the multinomial mixture model
  // from the mixture weights of the Poisson mixture model. 
  x %= u;
  x /= sum(x);

  // Perform one or more EM updates for the multinomial mixture model.
  mixem(L1,w,x,P,numiter);

  // Recover the mixture weights of the Poisson mixture model from the
  // mixture weights of the multinomial mixture model.
  x *= sum(w);
  x /= u;
}

// This third poismixem interface is similar to the second, except
// that the indices of the nonzero counts are supplied in vector i,
// and w must only contain the nonzero counts (vectors i and w should
// be the same length).
//
// For example, x1 and x2 should be the same after running this R code:
//
//   i  <- which(w > 0)
//   x1 <- poismixem2_rcpp(L1,w,u,x0,numiter)
//   x2 <- poismixem3_rcpp(L1,w[i],u,i-1,x0,numiter)
//
void poismixem (const mat& L1, const vec& u, const vec& w, const uvec& i,
		vec& x, uint numiter) {
  uint n = i.n_elem;
  uint k = x.n_elem;
  mat  P(n,k);
  poismixem(L1.rows(i),u,w,x,P,numiter);
}
