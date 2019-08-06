// This is included to suppress the warnings from solve() when the
// system is singular or close to singular.
#define ARMA_DONT_PRINT_ERRORS

#include <RcppArmadillo.h>
#include "misc.h"
#include "mixsqp.h"

using namespace Rcpp;
using namespace arma;

// FUNCTION DECLARATIONS
// ---------------------
void activesetqp(const mat& H, const vec& g, vec& y, int maxiter,
		 double zerothresholdsearchdir, double tol);
void compute_activeset_searchdir (const mat& H, const vec& y, vec& p, mat& B);
void backtracking_line_search (double f, const mat& L, const vec& w,
			       const vec& g, const vec& x, const vec& p,
			       const vec& e, double suffdecr, double beta,
			       double amin, vec& y, vec& u);
void   compute_grad (const mat& L, const vec& w, const vec& x,const vec& e,
		     vec& g, mat& H, mat& Z);
void   feasible_stepsize (const vec& x, const vec& p, int& j, double& a);
double init_hessian_correction (const mat& H, double a0);
double compute_objective (const mat& L, const vec& w, const vec& x,
			  const vec& e, vec& u);

// FUNCTION DEFINITIONS
// --------------------
// SQP algorithm for computing a maximum-likelihood estimate of a
// mixture model. For more information, see the help and comments
// accompanying the mixsqp R function and, in particular, see how
// mixsqp_rcpp is called inside the mixsqp function.
// 
// [[Rcpp::export]]
arma::vec mixsqp_rcpp (const arma::mat& L, const arma::vec& w,
		       const arma::vec& x0, double tol, double zerothreshold,
		       double zerosearchdir, double suffdecr,
		       double stepsizereduce, double minstepsize,
		       const arma::vec& e, int numiter, int maxiteractiveset,
		       bool verbose) {
  
  // Get the number of rows (n) and columns (m) of the conditional
  // likelihood matrix.
  int n = L.n_rows;
  int m = L.n_cols;

  // Initialize the solution.
  vec x = x0;
  
  // Scalars, vectors and matrices used in the computations below.
  double obj;
  uvec   i(m);
  vec    g(m);
  vec    ghat(m);
  vec    p(m);
  vec    u(n);
  mat    H(m,m);
  mat    Z(n,m);
  vec    y(m);
  
  // Iterate the SQP updates for a fixed number of iterations.
  for (int iter = 0; iter < numiter; iter++) {

    // Zero any co-ordinates that are below the specified threshold.
    i = find(x <= zerothreshold);
    x(i).fill(0);
    
    // Compute the value of the objective at x.
    obj = compute_objective(L,w,x,e,u);
    if (verbose)
      Rprintf("%4d %+0.15f\n",iter,obj);

    // Compute the gradient and Hessian.
    compute_grad(L,w,x,e,g,H,Z);
    
    // This is also a good point to check for a user interrupt; if the
    // user requests an interrupt, then an exception is thrown and
    // control is returned to the R console.
    Rcpp::checkUserInterrupt();
    
    // Solve the quadratic subproblem to obtain a search direction.
    ghat = g - H*x;
    y    = x;
    activesetqp(H,ghat,y,maxiteractiveset,zerosearchdir,tol);
    p = y - x;
    
    // Run backtracking line search.
    backtracking_line_search(obj,L,w,g,x,p,e,suffdecr,stepsizereduce,
			     minstepsize,y,u);
    x = y;
  }

  if (verbose) {
    obj = compute_objective(L,w,x,e,u);
    Rprintf("%4d %+0.15f\n",numiter,obj);
  }
  return x;
}

// This implements the active-set method from p. 472 of of Nocedal &
// Wright, Numerical Optimization, 2nd ed, 2006.
void activesetqp (const mat& H, const vec& g, vec& y, int maxiter,
		  double zerosearchdir, double tol) {
  int    m = g.n_elem;
  double a;
  int    k;
  vec    b(m);
  vec    p(m);
  vec    bs(m);
  vec    ps(m);
  vec    r(m);
  mat    Hs(m,m);
  mat    B(m,m);
  uvec   S(m);
  uvec   i(m);
  uvec   j(m);
  
  // This vector is used to keep track of the working set; all zero
  // entries of "t" are co-ordinates belonging to the working set.
  uvec t = (y > 0);
  
  // Run active set method to solve the quadratic subproblem.
  for (int iter = 0; iter < maxiter; iter++) {

    // Find the co-ordinates inside (j) and outside (i) the working
    // set.
    i = find(t != 0);
    j = find(t == 0);

    // Make sure that co-ordinates in the working set are set to zero.
    y(j).fill(0);

    // Define the smaller quadratic subproblem.
    Hs    = H(i,i);
    b     = g;
    b(i) += Hs*y(i);
    bs    = b(i);
      
    // Solve the quadratic subproblem to obtain a search direction.
    p.fill(0);
    compute_activeset_searchdir(Hs,bs,ps,B);
    p(i) = ps;
      
    // Check that the search direction is close to zero.
    if ((p.max() <= zerosearchdir) & (-p.min() <= zerosearchdir)) {
        
      // If all the Lagrange multiplers in the working set (that is,
      // zeroed co-ordinates) are positive, or nearly positive, we
      // have reached a suitable solution.
      if (j.is_empty())
	break;
      else if (b(j).min() >= -tol)
	break;
      else {
	
        // Find a co-ordinate with the smallest multiplier, and remove
        // it from the working set.
        k    = j(b(j).index_min());
        t(k) = 1;
      }
      
    // In this next part, we consider adding a co-ordinate to the
    // working set (but only if there are two or more non-zero
    // co-ordinates).
    } else {
        
      // Define the step size.
      feasible_stepsize(y,p,k,a);
      if (k >= 0 & a < 1) {
	
        // A blocking constraint exists; find it, and add it to the
        // working set (but only if there are two or more non-zero
        // co-ordinates).
	if (i.n_elem > 1)
	  t(k) = 0;
      }
      
      // Move to the new iterate along the search direction.
      y += a*p;
    }
  }
}

// Get the initial scalar multiplier for the identity matrix based on
// examining the diagonal entries of the Hessian.
double init_hessian_correction (const mat& H, double a0) {
  double d  = H.diag().min();
  double a;
  if (d > a0)
    a = 0;
  else
    a = a0 - d;
  return a;
}

// This implements Algorithm 3.3, "Cholesky with added multiple of the
// identity", from Nocedal & Wright, 2nd ed, p. 51.
void compute_activeset_searchdir (const mat& H, const vec& y, vec& p, mat& B) {
  double a0   = 1e-15;
  double amax = 1;
  double ainc = 10;
  int    n    = y.n_elem;
  mat    I(n,n,fill::eye);
  mat    R(n,n);

  // Get the initial scalar multiplier for the identity matrix.
  double a = init_hessian_correction(H,a0);

  // Repeat until a modified Hessian is found that is symmetric
  // positive definite, or until we cannot modify it any longer.
  while (true) {

    // Compute the modified Hessian.
    B = H + a*I;
    
    // Attempt to compute the Cholesky factorization of the modified
    // Hessian. If this fails, increase the contribution of the
    // identity matrix in the modified Hessian.
    if (chol(R,B))
      break;
    else if (a*ainc > amax)
      break;
    else if (a <= 0)
      a = a0;
    else
      a *= ainc;
  }
  
  // Compute the search direction using the modified Hessian.
  p = solve(B,-y);
}

// This implements the backtracking line search algorithm from p. 37
// of Nocedal & Wright, Numerical Optimization, 2nd ed, 2006.
void backtracking_line_search (double f, const mat& L, const vec& w,
			       const vec& g, const vec& x, const vec& p,
			       const vec& e, double suffdecr, double beta,
			       double amin, vec& y, vec& u) {
  int    i;
  double afeas;
  double fnew;

  // Determine the largest step size maintaining feasibility; if it is
  // larger than the minimum step size, return the minimum step size
  // that maintains feasibility of the solution. Otherwise, continue
  // to backtracking line search.
  feasible_stepsize(x,p,i,afeas);
  if (afeas <= amin)
    y = x + afeas*p;
  else {

    // Set the initial step size.
    double a = min(0.99,afeas);
    
    // Iteratively reduce the step size until either (1) we can't reduce
    // any more (because we have hit the minimum step size constraint),
    // or (2) the new candidate solution satisfies the "sufficient
    // decrease" condition.
    while (true) {
      y    = x + a*p;
      fnew = compute_objective(L,w,y,e,u);

      // Check whether the new candidate solution satisfies the
      // sufficient decrease condition, and remains feasible. If so,
      // accept this candidate solution.
      if (y.min() >= 0 & fnew <= f + suffdecr*a*dot(p,g))
        break;

      // If we cannot decrease the step size further, terminate the
      // backtracking line search, and set the step size to be the
      // minimum step size.
      else if (a*beta < amin) {
        y = x + amin*p;
        break;
      }
    
      // The new candidate does not satisfy the sufficient decrease
      // condition, so we need to try again with a smaller step size.
      a *= beta;
    }
  }
}

// Return the largest step size maintaining feasibility (x >= 0) for
// the given the search direction (p).
void feasible_stepsize (const vec& x, const vec& p, int& j, double& a) {
  uvec i = find(p < 0);
  a = 1;
  j = -1;
  if (!i.is_empty()) {
    vec t = -x(i)/p(i);
    j = t.index_min();
    if (t(j) < 1)
      a = t(j);
    j = i(j);
  }
}

// Compute the value of the objective at x; arguments L and w specify
// the objective, and e is a vector in which the entries can be set to
// small, positive numbers, or to zero. Input u stores an intermediate
// result used in the calculation.
double compute_objective (const mat& L, const vec& w, const vec& x,
			  const vec& e, vec& u) {
  u = L*x + e;
  if (u.min() <= 0)
    Rcpp::stop("Objective is -Inf");
  return sum(x) - sum(w % log(u));
}

// Compute the gradient and Hessian of the objective at x.
void compute_grad (const mat& L, const vec& w, const vec& x,
		  const vec& e, vec& g, mat& H, mat& Z) {
  vec u = L*x + e;
  g = -trans(L)*(w/u) + 1;
  Z = L;
  Z.each_col() %= (sqrt(w)/u);
  H = trans(Z) * Z;
}

