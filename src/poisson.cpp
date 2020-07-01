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
		     double& f1, vec& z0, vec& z1, vec& u, vec& loglik, 
		     double e, unsigned int numiter);

void fit_poisson_em_sparse (const vec& x, const vec& s, const vec& q,
			    double a, double b, double& f0, double& f1, 
			    vec& z0, vec& z1, vec& u, vec& loglik, double e, 
			    unsigned int numiter);

void fit_univar_poisson_models_em (const mat& X, const mat& L, const vec& s, 
				   mat& F0, mat& F1, mat& loglik, double e, 
				   unsigned int numiter, bool verbose);

// FUNCTION DEFINITIONS
// --------------------
// This implements fit_univar_poisson_models for method = "em-rcpp".
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
// [[Rcpp::export]]
List fit_univar_poisson_models_em_rcpp (const arma::mat& X, 
					const arma::mat& L,
					const arma::vec& s, double e,
					unsigned int numiter,
					bool verbose) {

  // Get the number of columns of the counts matrix (m) and the number
  // of topics (k).
  unsigned int m = X.n_cols;
  unsigned int k = L.n_cols;

  // Initialize the outputs.
  mat F0(m,k);
  mat F1(m,k);
  mat loglik(m,k);

  // Compute MLEs of the Poisson model parameters for each (j,k)
  // combination, where j is a row of the counts matrix, X, and k is a
  // topic (column of the L matrix).
  fit_univar_poisson_models_em(X,L,s,F0,F1,loglik,e,numiter,verbose);

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
  unsigned int n = x.n_elem;
  vec loglik(numiter);
  vec z0(n);
  vec z1(n);
  vec u(n);
  fit_poisson_em(x,s,q,f0,f1,z0,z1,u,loglik,e,numiter);
  return List::create(Named("f0")     = f0,
		      Named("f1")     = f1,
		      Named("loglik") = loglik);
}

// This is mainly used to test the fit_poisson_em function.
//
// [[Rcpp::export]]
List fit_poisson_em_sparse_rcpp (const arma::vec& x, const arma::vec& s,
				 const arma::vec& q, double a, double b, 
				 double f0, double f1, double e, 
				 unsigned int numiter) {
  unsigned int n = x.n_elem;
  vec loglik(numiter);
  vec z0(n);
  vec z1(n);
  vec u(n);
  fit_poisson_em_sparse(x,s,q,a,b,f0,f1,z0,z1,u,loglik,e,numiter);
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
// Vectors z0, z1 and u are used to store the E-step calculations.
// They should each be a vector of the same length as x.
//
// This should give the same result as the R function of the same name.
//
void fit_poisson_em (const vec& x, const vec& s, const vec& q, double& f0, 
		     double& f1, vec& z0, vec& z1, vec& u, vec& loglik, 
		     double e, unsigned int numiter) {

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

// This does the same thing as fit_poisson_em, except that the
// computations are better suited to spare matrices; see function
// fit_poisson_em for details. Input vector x need only contain the
// nonzero counts; inputs s and q should contain the entries
// corresponding to the nonzero counts in x. Input arguments z0, z1
// and u should be vectors of the same length as x. And input
// arguments a and b should be equal to a = sum(s*(1-q)) and b =
// sum(s*q), where s and q are the vectors corresponding to all
// counts, including the nonzero counts.
void fit_poisson_em_sparse (const vec& x, const vec& s, const vec& q,
			    double a, double b, double& f0, double& f1, 
			    vec& z0, vec& z1, vec& u, vec& loglik, double e, 
			    unsigned int numiter) {

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
    loglik[iter] = sum(x % log(u + e)) - (f0*a + f1*b);
  }
}

// Compute MLEs of the Poisson model parameters for each (j,k)
// combination, where j is a row of the counts matrix, X, and k is a
// topic (column of the L matrix). The parameter estimates are stored
// in F0, F1, and the log-likelihoods achieved at these estimates are
// stored in loglik; each of these should be m x k matrices, where m
// is the number of columns in the counts matrix, and k is the number
// of topics.
void fit_univar_poisson_models_em (const mat& X, const mat& L, const vec& s, 
				   mat& F0, mat& F1, mat& loglik, double e, 
				   unsigned int numiter, bool verbose) {

  // Get the number of rows (n) and columns (m) of the counts matrix, X,
  // and get the number of topics (k).
  unsigned int n = X.n_rows;
  unsigned int m = X.n_cols;
  unsigned int k = L.n_cols;

  // These variables are used to store the results from a single call
  // to fit_poisson_em.
  vec    ll(numiter);
  vec    z0(n);
  vec    z1(n);
  vec    u(n);
  double f0, f1;

  // Repeat for each column of the counts matrix, X, and for each
  // topic.
  Progress pb(m,verbose);
  for (unsigned int i = 0; i < m; i++) {
    pb.increment();
    checkUserInterrupt();
    for (unsigned int j = 0; j < k; j++) {
      f0 = 1;
      f1 = 1;
      fit_poisson_em(X.col(i),s,L.col(j),f0,f1,z0,z1,u,ll,e,numiter);
      F0(i,j)     = f0;
      F1(i,j)     = f1;
      loglik(i,j) = ll(numiter - 1);
    }
  }
}
