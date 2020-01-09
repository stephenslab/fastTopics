#include "poismix.h"
#include "misc.h"
#include "mixem.h"

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

// This is mainly used to test the first variant of the poismixsqp C++
// interface.
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec poismixsqp_rcpp (const arma::mat& L, const arma::vec& w,
			   const arma::vec& x0, uint numiter,
			   const Rcpp::List control) {
  mixsqp_control_params ctrl = get_mixsqp_control_params(control);
  return poismixsqp(L,w,x0,numiter,ctrl);
}

// This is mainly used to test the second variant of the poismixem C++
// interface.
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

// This is mainly used to test the second variant of the poismixsqp C++
// interface.
//
// [[Rcpp::export]]
arma::vec poismixsqp2_rcpp (const arma::mat& L1, const arma::vec& w,
			    const arma::vec& u, const arma::vec& x0,
			    uint numiter, const Rcpp::List control) {
  mixsqp_control_params ctrl = get_mixsqp_control_params(control);
  uint m = L1.n_cols;
  vec  x = x0;
  mat  Z = L1;
  mat  H(m,m);
  poismixsqp(L1,u,w,x,Z,H,numiter,ctrl);
  return x;
}

// This is mainly used to test the third variant of the poismixem C++
// interface.
//
// [[Rcpp::export]]
arma::vec poismixem3_rcpp (const arma::mat& L1, const arma::vec& w,
			   const arma::vec& u, const arma::uvec& i,
			   const arma::vec& x0, uint numiter) {
  vec x = x0;
  poismixem(L1,u,w,i,x,numiter);
  return x;
}

// This is mainly used to test the third variant of the poismixsqp C++
// interface.
//
// [[Rcpp::export]]
arma::vec poismixsqp3_rcpp (const arma::mat& L1, const arma::vec& w,
			    const arma::vec& u, const arma::uvec& i,
			    const arma::vec& x0, uint numiter,
			    const Rcpp::List control) {
  mixsqp_control_params ctrl = get_mixsqp_control_params(control);
  uint m = L1.n_cols;
  vec  x = x0;
  mat  H(m,m);
  poismixsqp(L1,u,w,i,x,H,numiter,ctrl);
  return x;
}

// Compute a maximum-likelihood estimate (MLE) of the mixture weights
// in a Poisson mixture model by iterating the multinomial mixture
// model EM updates for a fixed number of iterations.
//
// Input argument L is an n x m matrix with non-negative entries;
// input w is a vector of length n containing a count or "pseudocount"
// associated with each row of L; input argument x0 is the initial
// estimate of the mixture weights; and input "numiter" specifies the
// number of EM updates to perform.
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

// Compute a maximum-likelihood estimate (MLE) of the mixture weights
// in a Poisson mixture model by iterating the mix-SQP updates for a
// fixed number of iterations.
//
// Input argument L is an n x m matrix with non-negative entries;
// input w is a vector of length n containing a count or "pseudocount"
// associated with each row of L; input argument x0 is the initial
// estimate of the mixture weights; and input "numiter" specifies the
// number of dmix-SQPM updates to perform.
//
// The return value is a vector of length m containing the updated
// mixture weights.
vec poismixsqp (const mat& L, const vec& w, const vec& x0, uint numiter,
		const mixsqp_control_params& control) {
  int m  = L.n_cols;
  mat L1 = L;
  mat Z  = L;
  vec x  = x0;
  vec u  = sum(L,0);
  mat H(m,m);
  normalizecols(L1);
  poismixsqp(L1,u,w,x,Z,H,numiter,control);
  return x;
}

// Use this variant of poismixem if you plan on calling poismixem
// multiple times with the same matrix, L. In this call, input u should
// contain the column sums, u = colSums(L), in which L is the matrix
// prior to normalization, and matrix L1 is the normalized version of
// L in which each column sums to 1; that is, L1 = normalize.cols(L).
// Input P should be a matrix of the same size as L1.
//
// Note that in this variant of poismixem, L1 and w do not need to
// contain all the data; once u is correctly calculated, any rows of
// L1 associated with zero weights have no effect, so only the vector
// of nonzero weights w, and the rows of L1 associated with those
// weights, need to be supplied.
void poismixem (const mat& L1, const vec& u, const vec& w, vec& x, mat& P, 
		uint numiter) {
  double s = sum(w);
  
  // Recover the mixture proportions of the multinomial mixture model
  // from the mixture weights of the Poisson mixture model. 
  x %= u;
  x /= sum(x);

  // Perform one or more EM updates for the multinomial mixture model.
  mixem(L1,w,x,P,numiter);

  // Recover the mixture weights of the Poisson mixture model from the
  // mixture weights of the multinomial mixture model.
  x *= s;
  x /= u;
}

// Use this variant of poismixsqp if you plan on calling poismixsqp
// multiple times with the same matrix, L. In this call, input u should
// contain the column sums, u = colSums(L), in which L is the matrix
// prior to normalization, and matrix L1 is the normalized version of
// L in which each column sums to 1; that is, L1 = normalize.cols(L).
// Input Z should be a matrix of the same size as L1, and H should be
// an m x m matrix, where m is the number of columns in L.
//
// Note that in this variant of poismixem, L1 and w do not need to
// contain all the data; once u is correctly calculated, any rows of
// L1 associated with zero weights have no effect, so only the vector
// of nonzero weights w, and the rows of L1 associated with those
// weights, need to be supplied.
void poismixsqp (const mat& L1, const vec& u, const vec& w, vec& x, mat& Z,
		 mat& H, uint numiter, const mixsqp_control_params& control) {
  uvec i = find(w > 0);
  if (i.n_elem == 1) {

    // Treat the special case when only one of the counts is nonzero.
    poismix_one_nonzero(L1,u,w,i(0),x);
  } else {
    vec    objective(numiter);
    double s = sum(w);

    // Recover the mixture proportions of the multinomial mixture model
    // from the mixture weights of the Poisson mixture model.
    x %= u;
    x /= s;

    // Perform one or more SQP updates for the multinomial mixture model.
    mixsqp(L1,w,x,Z,H,numiter,control,objective);

    // Recover the mixture weights of the Poisson mixture model from the
    // mixture weights of the multinomial mixture model.
    x *= s;
    x /= u;
  }
}

// This third poismixem interface is similar to the second, except
// that the indices of the nonzero counts are supplied in vector i,
// and w must contain only the nonzero counts (vectors i and w should
// be the same length).
//
// To illustrate, in this example x1 and x2 should be the same after
// running this R code:
//
//   i  <- which(w > 0)
//   x1 <- poismixem2_rcpp(L1,w,u,x0,numiter)
//   x2 <- poismixem3_rcpp(L1,w[i],u,i-1,x0,numiter)
//
void poismixem (const mat& L1, const vec& u, const vec& w, const uvec& i,
		vec& x, uint numiter) {
  uint n = i.n_elem;
  uint m = x.n_elem;
  mat  P(n,m);
  poismixem(L1.rows(i),u,w,x,P,numiter);
}

// The third poismixsqp interface is similar to the second, except
// that the indices of the nonzero counts are supplied in vector i,
// and w must contain only the nonzero counts (i and w should be the
// same length).
//
// To illustrate, in this example x1 and x2 should be the same after
// running this R code:
//
//   i  <- which(w > 0)
//   x1 <- poismixsqp2_rcpp(L1,w,u,x0,numiter,control)
//   x2 <- poismixsqp3_rcpp(L1,w[i],u,i-1,x0,numiter,control)
//
void poismixsqp (const mat& L1, const vec& u, const vec& w, const uvec& i,
		 vec& x, mat& H, uint numiter,
		 const mixsqp_control_params& control) {
  uint n = i.n_elem;
  uint m = x.n_elem;
  mat  Z(n,m);
  vec  objective(numiter);
  poismixsqp(L1.rows(i),u,w,x,Z,H,numiter,control);
}

// Find the maximum-likelihood estimate (MLE) for the special case
// when only one of the counts is positive. Input argument w should be
// a vector of counts (in which exactly one of these counts is
// positive), and i should be the index of the nonzero count. See
// above for an explanation of inputs L1 and u.
void poismix_one_nonzero (const mat& L1, const vec& u, const vec& w, 
			  uint i, vec& x) {
  mixture_one_nonzero(L1,i,x);
  uint j = max(x);
  x(j)   = w(i)/u(j);
}
