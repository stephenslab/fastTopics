#include "misc.h"

using namespace Rcpp;
using namespace arma;

// FUNCTION DECLARATIONS
// ---------------------
double loglik_poisson (const vec& x, const vec& y, double sumy, double e);

double simulate_posterior_poisson (const vec& x, const mat& L, const vec& w,
				   const vec& f0, const mat& D, const mat& U, 
				   mat& samples, double s, double e);

// FUNCTION DEFINITIONS
// --------------------
// This is a more efficient C++ implementation of
// simulate_posterior_poisson. See the description of that function
// for an explanation of the input arguments. Additional input
// arguments specific to the C++ implementation are D and U; these
// inputs should be generated in R as
//
//   k = length(f) 
//   D = matrix(rnorm(ns*k),ns,k)
//   U = matrix(runif(ns*k),ns,k)
//
// where ns is the requesed number of Monte Carlo samples.
// 
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List simulate_posterior_poisson_rcpp (const arma::vec& x, const arma::mat& L,
				      const arma::vec& f, const arma::mat& D,
				      const arma::mat& U, double s, double e) {
  unsigned int ns = D.n_rows;
  unsigned int k  = D.n_cols;
  mat samples(ns,k);
  double ar = simulate_posterior_poisson(x,L,sum(L,0),f,D,U,samples,s,e);
  return List::create(Named("samples") = samples,Named("ar") = ar);
}

// This is the same as simulate_posterior_poisson_rcpp, except that w
// = colSums(L) is precomputed. This will produce the correct result
// when the zero counts x are omitted, which makes this particularly
// useful when the counts are sparse.
//
// [[Rcpp::export]]
List simulate_posterior_poisson_sparse_rcpp (const arma::vec& x, 
					     const arma::mat& L,
					     const arma::vec& w, 
					     const arma::vec& f,
					     const arma::mat& D,
					     const arma::mat& U, 
					     double s, double e) {
  unsigned int ns = D.n_rows;
  unsigned int k  = D.n_cols;
  mat samples(ns,k);
  double ar = simulate_posterior_poisson(x,L,w,f,D,U,samples,s,e);
  return List::create(Named("samples") = samples,Named("ar") = ar);
}

// This should give the same, or nearly the same, result as
// sum(dpois(x,y,log = TRUE)) so long as sumy = sum(y), except that
// terms that do not depend on the Poisson rates u are not included.
double loglik_poisson (const vec& x, const vec& y, double sumy, double e) {
  return sum(x % log(y + e)) - sumy;
}

// This implements the core part of simulate_posterior_poisson_rcpp
// and simulate_posterior_poisson_sparse_rcpp. Input argument
// "samples" should be an ns x k matrix, where ns is the number of
// Monte Carlo samples to generate, and k = length(t).
double simulate_posterior_poisson (const vec& x, const mat& L, const vec& w,
				   const vec& f0, const mat& D, const mat& U, 
				   mat& samples, double s, double e) {
  unsigned int n  = x.n_elem;
  unsigned int ns = samples.n_rows;
  unsigned int k  = samples.n_cols;
  vec t = log(f0);
  vec f(k);
  vec fnew(k);
  vec tnew(k);
  vec u(n);
  vec unew(n);
  double ar = 0;
  double ll, llnew;
  double a, d;
  for (unsigned int i = 0; i < ns; i++) {
    u = L * exp(t);
    for (unsigned int j = 0; j < k; j++) {

      // Randomly suggest moving to tj(new) = tj + d, where d ~ N(0,s).
      d        = s*D(i,j);
      tnew     = t;
      tnew(j) += d;

      // Compute the Metropolis acceptance probability, and move to the
      // new state according to this acceptance probability. Note that
      // the additional d in the acceptance probability is needed to
      // account for the fact that we are simulating log(f), not f; see
      // p. 11 of Devroye (1986) "Non-uniform random variate generation".
      f     = exp(t);
      fnew  = exp(tnew);
      unew  = u + (exp(tnew[j]) - exp(t[j])) * L.col(j);
      ll    = loglik_poisson(x,u,sum(w % f),e);
      llnew = loglik_poisson(x,unew,sum(w % fnew),e);
      a     = exp((llnew - ll) + d);
      a     = minimum(1,a);
      if (U(i,j) < a) {
        t = tnew;
	u = unew;
        ar++;
      }
    }

    // Store the current state of the Markov chain.
    samples.row(i) = exp(t);
  }

  ar /= k*ns;
  return ar;
}

