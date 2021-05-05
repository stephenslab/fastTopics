#include "misc.h"

using namespace Rcpp;
using namespace arma;

// FUNCTION DECLARATIONS
// ---------------------
double loglik_poisson (const vec& x, const vec& y, double e);

double simulate_posterior_poisson (const vec& x, const mat& L, const mat& D, 
				   const mat& U, vec& t, mat& samples, 
				   double s, double e);

// FUNCTION DEFINITIONS
// --------------------
// TO DO: Explain here what this function does, and how to use it.
// 
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List simulate_posterior_poisson_rcpp (const arma::vec& x, const arma::mat& L,
				      const arma::vec& f, const arma::mat& D,
				      const arma::mat& U, double s, double e) {
  unsigned int ns = D.n_rows;
  unsigned int k  = D.n_cols;
  vec t = log(f);
  mat samples(ns,k);
  double ar = simulate_posterior_poisson(x,L,D,U,t,samples,s,e);
  return List::create(Named("samples") = samples,Named("ar") = ar);
}

// This should give the same, or nearly the same, result as
// sum(dpois(x,y,log = TRUE)), except that terms that do not depend on
// the Poisson rates u are not included
double loglik_poisson (const vec& x, const vec& y, double e) {
  return sum(x % log(y + e) - y);
}

// TO DO: Explain here what this function does, and how to use it.
double simulate_posterior_poisson (const vec& x, const mat& L, const mat& D, 
				   const mat& U, vec& t, mat& samples, 
				   double s, double e) {
  unsigned int n  = x.n_elem;
  unsigned int ns = samples.n_rows;
  unsigned int k  = samples.n_cols;
  vec    tnew(k);
  vec    u(n);
  vec    unew(n);
  double ar = 0;
  double ll, llnew;
  double a, d;
  for (unsigned int i = 0; i < ns; i++) {
    u = L * exp(t);
    for (unsigned int j = 0; j < k; j++) {

      // Randomly suggest moving to tj(new) = tj + d, where d ~ N(0,s).
      tnew = t;
      d = s*D(i,j);
      tnew(j) += d;

      // Compute the Metropolis acceptance probability, and move to the
      // new state according to this acceptance probability. Note that
      // the additional d in the acceptance probability is needed to
      // account for the fact that we are simulating log(f), not f; see
      // p. 11 of Devroye (1986) "Non-uniform random variate generation".
      unew  = u + (exp(tnew[j]) - exp(t[j])) * L.col(j);
      ll    = loglik_poisson(x,u,e);
      llnew = loglik_poisson(x,unew,e);
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

