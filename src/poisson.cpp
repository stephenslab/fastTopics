#include <RcppArmadillo.h>
#include <progress.hpp>
#include <progress_bar.hpp>

using namespace Rcpp;
using namespace arma;

// FUNCTION DECLARATIONS
// ---------------------
void get_poisson_rates (const vec& s, const vec& q, double f0, double f1, 
			vec& su);

double loglik_poisson (const vec& x, const vec& u, double e);

void fit_poisson_em (const vec& x, const vec& s, const vec& q, double& f0, 
		     double& f1, vec& loglik, double e, unsigned int numiter);

// FUNCTION DEFINITIONS
// --------------------
// TO DO: Explain here what this function does, and how to use it.
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
// [[Rcpp::export]]
List fit_univar_poisson_models_em_rcpp (const arma::mat& X, 
					const arma::mat& L,
					const arma::vec& s, double e,
					unsigned int numiter,
					bool verbose) {

  // Get the number of columns of the counts matrix (m), and the
  // number of topics (k).
  unsigned int m = X.n_cols;
  unsigned int k = L.n_cols;

  // Initialize the outputs.
  mat F0(m,k,fill::ones);
  mat F1(m,k,fill::ones);
  mat loglik(m,k);

  // Repeat for each column of the counts matrix, X, and for each
  // topic.
  double f0, f1;
  for (unsigned int i = 0; i < m; i++) {
    for (unsigned int j = 0; j < k; j++) {
      // loglik(i,j) = fit_poisson_em(X[,i],s,L[,j],e = e,numiter = numiter)
    }
  }
  
  // Output the MLEs of the Poisson model parameters (F0, F1), and the
  // values of the Poisson model log-likelihood attained at those
  // parameters.
  return List::create(Named("F0")     = F0,
                      Named("F1")     = F1,
		      Named("loglik") = loglik);
}

// This is mainly used to test the fit_poisson_em function.
//
// [[Rcpp::export]]
List fit_poisson_em_rcpp (const arma::vec& x, const arma::vec& s,
			  const arma::vec& q, double f0, double f1,
			  double e, unsigned int numiter) {
  vec loglik(numiter);
  fit_poisson_em(x,s,q,f0,f1,loglik,e,numiter);
  return List::create(Named("f0")     = f0,
		      Named("f1")     = f1,
		      Named("loglik") = loglik);
}

// Compute the Poisson rates given the size factors (s), topic
// proportions (q), and parameters (f0, f1) of the Poisson model.
void get_poisson_rates (const vec& s, const vec& q, double f0, double f1, 
			vec& su) {
  su = s % (f0*(1-q) + f1*q);
}

// This should give the same, or nearly the same, result as
// sum(dpois(x,u,log = TRUE)), except that terms that do not depend on
// the Poisson rates (u) are disregarded.
double loglik_poisson (const vec& x, const vec& u, double e) {
  double y = sum(x % log(u + e)) - sum(u);
  return y;
}

// Perform EM updates to compute maximum-likelihood estimates (MLEs)
// of the parameters in the single-count Poisson model: x ~
// Poisson(s*u), with u given by u = f0*(1-q) + f1*q. Parameters f0,
// f1 are estimated, and vectors s, q are provided. Input arguments f0
// and f1 are initial estimates of the parameters. Input argument "e"
// is a small positive scalar added to the denominator in the E-step
// to avoid division by zero.
//
// This should give the same result as the R function of the same name.
//
void fit_poisson_em (const vec& x, const vec& s, const vec& q, double& f0, 
		     double& f1, vec& loglik, double e, unsigned int numiter) {

  // Get the number of samples.
  unsigned int n = x.n_elem;

  // These vectors are used to store the posterior expectations
  // computed in the E-step.
  vec z0(n);
  vec z1(n);
  vec u(n);

  // Store a couple pre-calculations to simplify the calculations
  // below.
  double a = sum(s % (1-q));
  double b = sum(s % q);

  // Perform the EM updates. Progress is monitored by computing the
  // log-likelihood at each iteration.
  for (unsigned int iter = 0; iter < numiter; iter++) {

    // E-STEP
    // ------
    z0  = f0*(1-q);
    z1  = f1*q;
    u   = z0 + z1 + e;
    z0 %= x/u;
    z1 %= x/u;

    // M-STEP
    // ------
    f0 = sum(z0)/a;
    f1 = sum(z1)/b;

    // Compute the log-likelihood at the current estimates of the
    // model parameters (ignoring terms that don't depend on f0 or f1).
    get_poisson_rates(s,q,f0,f1,u);
    loglik[iter] = loglik_poisson(x,u,e);
  }
}
