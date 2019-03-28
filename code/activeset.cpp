// This is included to suppress the warnings from solve() when the
// system is singular or close to singular.
#define ARMA_DONT_PRINT_ERRORS

#include <RcppArmadillo.h>

// This depends statement is needed to tell R where to find the
// additional header files.
//
// [[Rcpp::depends(RcppArmadillo)]]
//

using namespace Rcpp;
using namespace arma;

// TO DO: Explain here what this function does, and how to use it.
//
// [[Rcpp::export]]
//
double activeset_rcpp (const mat& H, const vec& g, vec& y,
		       int maxiter_activeset, double zerothreshold) {

  // Get the number of parameters to optimize.
  int k = y.n_elem;

  double n;
  uvec   t(k);
  uvec   i(k);
  vec    b(k);
  vec    bs(k);
  mat    Hs(k,k);

  // Initialize the solution to the quadratic subproblem.
  t = (y >= zerothreshold);
  n = sum(t);
  y.fill(0);
  y.elem(i).fill(1/n);
  
  // Run active set method to solve the quadratic subproblem.
  for (int iter = 0; iter < maxiter_activeset; iter++) {

    // Define the equality-constrained quadratic subproblem.
    b  = H*y + g;
    i  = find(t);
    bs = b.elem(i);
    Hs = H.elem(i,i);
  }
  
  return g.min();
}
