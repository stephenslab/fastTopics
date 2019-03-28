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
vec activeset_rcpp (const mat& H, const vec& g, const vec& y0,
		    int maxiter_activeset, double convtol,
		    double zerothreshold, double zerosearchdir) {

  // Get the number of parameters to be optimized.
  int k = g.n_elem;

  double n;
  double a;
  int    inew;
  
  vec  y = y0;
  uvec t(k);
  uvec i(k);
  uvec j(k);
  uvec A(k);
  vec  p(k);
  vec  p0(k);
  vec  alpha(k);
  vec  b(k);
  vec  bs(k);
  mat  Hs(k,k);

  // Initialize the solution to the quadratic subproblem.
  t = (y >= zerothreshold);
  n = sum(t);
  y.fill(0);
  y.elem(i).fill(1/n);
  
  // Run active set method to solve the quadratic subproblem.
  for (int iter = 0; iter < maxiter_activeset; iter++) {

    // Get the set of co-ordinates outside the working set.
    i = find(t);
    j = find(1 - t);
    n = j.n_elem;
    
    // Define the equality-constrained quadratic subproblem.
    b  = H*y + g;
    bs = b.elem(i);
    Hs = H.elem(i,i);

    // Solve the equality-constrained subproblem.
    p.fill(0);
    p.elem(i) = -solve(Hs,bs);
      
    // Reset the step size.
    a = 0.99;

    // If the working set is empty, and we have already tried to
    // update the working set at least once, we have reached a
    // suitable solution.
    if (n == 0 & iter > 0) {
      break;

    // Check that the search direction is close to zero.
    } else if ((p.max()  <= zerosearchdir) &
	       (-p.min() <= zerosearchdir) &
	       (n > 0)) {

      // If all the gradient entries in the working set (that is,
      // zeroed co-ordinates) are positive, or nearly positive, we
      // have reached a suitable solution.
      if (b(j).min() >= convtol)
	break;

      // Find a co-ordinate with the smallest gradient entry, and
      // remove it from the working set.
      inew    = j[b(j).index_min()];
      t[inew] = 1;

     // In this next part, we consider adding a co-ordinate to the
     // working set, but only if there are two or more non-zero
     // co-ordinates.
    } else if (n < k - 1) {

      // Revise the step size.
      p0 = p;
      p0.elem(j).fill(0);
      A = find(p0 < 0);
      if (!A.is_empty()) {
        alpha = -y.elem(A)/p.elem(A);
        inew  = alpha.index_min();
        if (alpha[inew] < 1) {

	  // Blocking constraint exists; find and add it to the
          // working set.
          a          = alpha[inew]; 
          t[A[inew]] = 0;
        }
      }
    }

    // Move to the new iterate along the search
    // direction.
    y += a * p;
  }
  
  return y;
}
