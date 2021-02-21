#include <RcppArmadillo.h>
#include <progress.hpp>
#include <progress_bar.hpp>
#include "misc.h"

using namespace Rcpp;
using namespace arma;

// FUNCTION DECLARATIONS
// ---------------------
void get_poisson_rates (const vec& s, const vec& q, double f0, double f1, 
			vec& su);

double loglik_poisson (const vec& x, const vec& u, double e);

double fit_poisson_em (const vec& x, const vec& s, const vec& q, double& f0, 
		       double& f1, vec& z0, vec& z1, vec& u, double e,
		       unsigned int& numiter, double tol);

double fit_poisson_em_sparse (const vec& x, const vec& s, const vec& q,
			      double a, double b, double& f0, double& f1, 
			      vec& z0, vec& z1, vec& u, double e,
			      unsigned int& numiter, double tol);

void fit_univar_poisson_models_em (const mat& X, const mat& L, const vec& s, 
				   mat& F0, mat& F1, mat& loglik, double e, 
				   unsigned int numiter, double tol,
				   bool verbose);

void fit_univar_poisson_models_em_sparse (const sp_mat& X, const mat& L, 
					  const vec& s, mat& F0, mat& F1, 
					  mat& loglik, double e, 
					  unsigned int numiter, double tol,
					  bool verbose);

// FUNCTION DEFINITIONS
// --------------------
// This implements fit_univar_poisson_models for method = "em-rcpp",
// with *dense* counts matrix X.
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
// [[Rcpp::export]]
List fit_univar_poisson_models_em_rcpp (const arma::mat& X, 
					const arma::mat& L,
					const arma::vec& s, double e,
					unsigned int numiter, double tol,
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
  fit_univar_poisson_models_em(X,L,s,F0,F1,loglik,e,numiter,tol,verbose);

  // Output the MLEs of the Poisson model parameters (F0, F1), and the
  // values of the Poisson model log-likelihood attained at those
  // parameters.
  return List::create(Named("F0")      = F0,
                      Named("F1")      = F1,
		      Named("loglik")  = loglik,
		      Named("numiter") = numiter);
}

// This implements fit_univar_poisson_models for method = "em-rcpp",
// with *sparse* counts matrix X.
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
// [[Rcpp::export]]
List fit_univar_poisson_models_em_sparse_rcpp (const arma::sp_mat& X, 
					       const arma::mat& L,
					       const arma::vec& s, double e,
					       unsigned int numiter,
					       double tol, bool verbose) {

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
  fit_univar_poisson_models_em_sparse(X,L,s,F0,F1,loglik,e,numiter,tol,
				      verbose);

  // Output the MLEs of the Poisson model parameters (F0, F1), and the
  // values of the Poisson model log-likelihood attained at those
  // parameters.
  return List::create(Named("F0")      = F0,
                      Named("F1")      = F1,
		      Named("loglik")  = loglik,
		      Named("numiter") = numiter);
}

// This is mainly used to test the fit_poisson_em function.
//
// [[Rcpp::export]]
List fit_poisson_em_rcpp (const arma::vec& x, const arma::vec& s,
			  const arma::vec& q, double f0, double f1,
			  double e, unsigned int numiter, double tol) {
  unsigned int n = x.n_elem;
  vec z0(n);
  vec z1(n);
  vec u(n);
  unsigned int t = numiter;
  double loglik = fit_poisson_em(x,s,q,f0,f1,z0,z1,u,e,t,tol);
  return List::create(Named("f0")      = f0,
		      Named("f1")      = f1,
		      Named("loglik")  = loglik,
		      Named("numiter") = t);
}

// This is mainly used to test the fit_poisson_em function.
//
// [[Rcpp::export]]
List fit_poisson_em_sparse_rcpp (const arma::vec& x, const arma::vec& s,
				 const arma::vec& q, double a, double b, 
				 double f0, double f1, double e, 
				 unsigned int numiter, double tol) {
  unsigned int n = x.n_elem;
  vec z0(n);
  vec z1(n);
  vec u(n);
  unsigned int t = numiter;
  double loglik = fit_poisson_em_sparse(x,s,q,a,b,f0,f1,z0,z1,u,e,t,tol);
  return List::create(Named("f0")      = f0,
		      Named("f1")      = f1,
		      Named("loglik")  = loglik,
		      Named("numiter") = t);
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
// The return value is the log-likelihood at the current estimates,
// ignoring terms that don't depend on f0 or f1. Also, the "numiter"
// input is updated to reflect the number of EM updates that were
// actually performed.
double fit_poisson_em (const vec& x, const vec& s, const vec& q, double& f0, 
		       double& f1, vec& z0, vec& z1, vec& u, double e,
		       unsigned int& numiter, double tol) {
  
  // Store a couple pre-calculations to simplify the calculations
  // below.
  double a = sum(s % (1-q)) + e;
  double b = sum(s % q) + e;
  double f00, f10;
  
  // Perform the EM updates. Progress is monitored by computing the
  // log-likelihood at each iteration.
  for (unsigned int iter = 0; iter < numiter; iter++) {

    // Save the current parameter estimates.
    f00 = f0;
    f10 = f1;
    
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

    // Stop early if the parameters have not changed much.
    if ((f0 - f00) < tol && 
	(f00 - f0) < tol && 
	(f1 - f10) < tol &&
	(f10 - f1) < tol) {

      // Update "numiter" to reflect the number of EM updates actually
      // performed.
      numiter = iter + 1;
      break;
    }
  }
  
  // Compute the log-likelihood, ignoring terms that don't depend on
  // f0 or f1.
  get_poisson_rates(s,q,f0,f1,u);
  return loglik_poisson(x,u,e);
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
//
// The return value is the log-likelihood at the current estimates,
// ignoring terms that don't depend on f0 or f1. Also, the "numiter"
// input is updated to reflect the number of EM updates that were
// actually performed.
double fit_poisson_em_sparse (const vec& x, const vec& s, const vec& q,
			      double a, double b, double& f0, double& f1, 
			      vec& z0, vec& z1, vec& u, double e, 
			      unsigned int& numiter, double tol) {
  double f00, f10;

  // Perform the EM updates. Progress is monitored by computing the
  // log-likelihood at each iteration.
  for (unsigned int iter = 0; iter < numiter; iter++) {

    // Save the current parameter estimates.
    f00 = f0;
    f10 = f1;
    
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

    // Stop early if the parameters have not changed much.
    if ((f0 - f00) < tol && 
	(f00 - f0) < tol && 
	(f1 - f10) < tol &&
	(f10 - f1) < tol) {

      // Update "numiter" to reflect the number of EM updates actually
      // performed.
      numiter = iter + 1;
      break;
    }
  }

  // Compute the log-likelihood, ignoring terms that don't depend on
  // f0 or f1.
  get_poisson_rates(s,q,f0,f1,u);
  return sum(x % log(u + e)) - (f0*a + f1*b);
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
				   unsigned int numiter, double tol,
				   bool verbose) {

  // Get the number of rows (n) and columns (m) of the counts matrix, X,
  // and get the number of topics (k).
  unsigned int n = X.n_rows;
  unsigned int m = X.n_cols;
  unsigned int k = L.n_cols;

  // These variables are used to store the results from a single call
  // to fit_poisson_em.
  vec    z0(n);
  vec    z1(n);
  vec    u(n);
  double f0, f1;
  unsigned int t;
  
  // Repeat for each column of the counts matrix, X, and for each
  // topic.
  Progress pb(m,verbose);
  for (unsigned int i = 0; i < m; i++) {
    pb.increment();
    checkUserInterrupt();
    for (unsigned int j = 0; j < k; j++) {
      f0 = 0.5;
      f1 = 0.5;
      t  = numiter;
      loglik(i,j) = fit_poisson_em(X.col(i),s,L.col(j),f0,f1,z0,z1,u,e,t,tol);
      F0(i,j) = f0;
      F1(i,j) = f1;
    }
  }
}

// This implements the EM algorithm for the case when the counts
// matrix, X, is sparse.
void fit_univar_poisson_models_em_sparse (const sp_mat& X, const mat& L, 
					  const vec& s, mat& F0, mat& F1, 
					  mat& loglik, double e, 
					  unsigned int numiter, double tol,
					  bool verbose) {

  // Get the number of rows (n) and columns (m) of the counts matrix, X,
  // and get the number of topics (k).
  unsigned int n = X.n_rows;
  unsigned int m = X.n_cols;
  unsigned int k = L.n_cols;

  // Pre-calculate, for each topic (j), the "a" and "b" inputs
  // required for fit_poisson_em_sparse.
  vec a(k);
  vec b(k);
  vec q(n);
  for (unsigned int j = 0; j < k; j++) {
    q    = L.col(j);
    a(j) = sum(s % (1-q)) + e;
    b(j) = sum(s % q) + e;
  }

  // Repeat for each column of the counts matrix, X, and for each
  // topic.
  Progress pb(m,verbose);
  double f0, f1;
  unsigned int t;
  for (unsigned int i = 0; i < m; i++) {

    // Extract the count data from the ith column, and set up other
    // data structures for running fit_poisson_em_sparse.
    vec          x  = nonzeros(X.col(i));
    unsigned int n1 = x.n_elem;
    uvec r(n1);
    getcolnonzeros(X,r,i);
    vec  si = s(r);
    vec  qi(n1);
    vec  z0(n1);
    vec  z1(n1);
    vec  u(n1);
    pb.increment();
    checkUserInterrupt();
    for (unsigned int j = 0; j < k; j++) {

      // Extract the topic proportions corresponding to the nonzero counts.
      getcolelems(L,r,j,qi);

      // Perform the EM updates.
      f0 = 0.5;
      f1 = 0.5;
      t  = numiter;
      loglik(i,j) = fit_poisson_em_sparse(x,si,qi,a(j),b(j),f0,f1,z0,z1,u,e,t,
					  tol);
      F0(i,j) = f0;
      F1(i,j) = f1;
    }
  }
}
