#include "misc.h"

using namespace Rcpp;
using namespace arma;

// FUNCTION DECLARATIONS
// ---------------------
double loglik_poisson (const vec& x, const vec& y, double sumy, double e);

double simulate_posterior_poisson (const vec& x, const mat& L, const vec& w,
				   const vec& f, const mat& D, const mat& U, 
				   const mat& M, mat& samples, double s, 
				   double e);

// FUNCTION DEFINITIONS
// --------------------
// This is a more efficient C++ implementation of the R function
// simulate_posterior_poisson. See the description of that function
// for an explanation of the input arguments. Additional input
// arguments specific to the C++ implementation are D, U and M; these
// inputs should be generated in R as
//
//   k = length(f) 
//   D = matrix(rnorm(ns*k),ns,k)
//   U = matrix(runif(ns*k),ns,k)
//   M = matrix(sample(k,ns*k,replace = TRUE),ns,k) - 1
//
// where ns is the number of Monte Carlo samples to simulate. Note in
// particular that the 1 needs to be subtracted from the random
// indices M because the indices in C++ start at zero.
// 
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List simulate_posterior_poisson_rcpp (const arma::vec& x, const arma::mat& L,
				      const arma::vec& f, const arma::mat& D,
				      const arma::mat& U, const arma::mat& M,
				      double s, double e) {
  unsigned int ns = D.n_rows;
  unsigned int k  = D.n_cols;
  mat samples(ns,k);
  double ar = simulate_posterior_poisson(x,L,sum(L,0),f,D,U,M,samples,s,e);
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
					     const arma::mat& M,
					     double s, double e) {
  unsigned int ns = D.n_rows;
  unsigned int k  = D.n_cols;
  mat samples(ns,k);
  double ar = simulate_posterior_poisson(x,L,w,f,D,U,M,samples,s,e);
  return List::create(Named("samples") = samples,Named("ar") = ar);
}

// This implements the core part of simulate_posterior_poisson_rcpp
// and simulate_posterior_poisson_sparse_rcpp. Input argument
// "samples" should be an ns x k matrix, where ns is the number of
// Monte Carlo samples to generate, and k = length(f).
double simulate_posterior_poisson (const vec& x, const mat& L, const vec& w,
				   const vec& f, const mat& D, const mat& U, 
				   const mat& M, mat& samples, double s, 
				   double e) {
  unsigned int n  = x.n_elem;
  unsigned int ns = samples.n_rows;
  unsigned int k  = samples.n_cols;
  unsigned int j;
  vec g = log(f);
  vec gnew(k);
  vec u(n);
  vec unew(n);
  double ar = 0;
  double ll, llnew, su, sunew, a, b, d;

  // Compute the log-density at the initial state.
  u  = L*exp(g);
  su = dot(w,exp(g));
  ll = loglik_poisson(x,u,su,e);
  for (unsigned int i = 0; i < ns; i++) {

    // Recompute the Poisson rates (u) and their sum (su).
    u  = L*exp(g);
    su = dot(w,exp(g));
    for (unsigned int t = 0; t < k; t++) {

      // Randomly suggest moving to gj* = gj + d, where d ~ N(0,s).
      j        = M(i,t);
      d        = s*D(i,t);
      gnew     = g;
      gnew(j) += d;

      // Compute the Metropolis acceptance probability, and move to
      // the new state according to this acceptance probability. Note
      // that the additional "d" in the acceptance probability is
      // needed to account for the fact that we are simulating g =
      // log(f), not f; see p. 11 of Devroye (1986) "Non-uniform
      // random variate generation".
      b     = exp(gnew[j]) - exp(g[j]);
      unew  = u + b*L.col(j);
      sunew = su + b*w[j];
      llnew = loglik_poisson(x,unew,sunew,e);
      a = exp((llnew - ll) + d);
      a = minimum(1,a);
      if (U(i,t) < a) {
        g  = gnew;
	u  = unew;
	su = sunew;
	ll = llnew;
        ar++;
      }
    }

    // Store the current state of the Markov chain.
    samples.row(i) = g;
  }

  ar /= k*ns;
  return ar;
}

// This should produce the same, or nearly the same, result as
// sum(dpois(x,y,log = TRUE)) so long as sumy = sum(y), except that
// terms that do not depend on the Poisson rates u are not included.
double loglik_poisson (const vec& x, const vec& y, double sumy, double e) {
  return dot(x,log(y + e)) - sumy;
}
