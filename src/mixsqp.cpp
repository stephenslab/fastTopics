#include "misc.h"
#include "mixsqp.h"

using namespace Rcpp;
using namespace arma;

// FUNCTION DECLARATIONS
// ---------------------
double compute_objective (const mat& L, const vec& w, const vec& x, double e);
void   compute_grad      (const mat& L, const vec& w, const vec& x, double e,
			  vec& g, mat& H, mat& Z);
void   activesetqp       (const mat& H, const vec& g, vec& y, int maxiter,
			  double zerosearchdir, double tol, double ainc);
void   compute_activeset_searchdir (const mat& H, const vec& y, vec& p, mat& B,
				    double ainc);
void   backtracking_line_search (double f, const mat& L, const vec& w,
				 const vec& g, const vec& x, const vec& y,
				 double e, double suffdecr, double beta,
				 double amin, vec& xnew);

// FUNCTION DEFINITIONS
// --------------------
// This is mainly used for testing the mixsqp C++ function.
// 
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List mixsqp_rcpp (const arma::mat& L, const arma::vec& w,
			const arma::vec& x0, unsigned int numiter,
			const Rcpp::List& control) {
  mixsqp_control_params ctrl = get_mixsqp_control_params(control);
  vec objective(numiter);
  vec x = mixsqp(L,w,x0,numiter,ctrl,objective);
  return List::create(Named("x") = x,Named("objective") = objective);
}

// Get the mix-SQP optimization settings from a named list in R.
mixsqp_control_params get_mixsqp_control_params	(const Rcpp::List& control) {
  mixsqp_control_params out;
  out.convtolactiveset        = control["convtol.activeset"];
  out.zerothresholdsolution   = control["zero.threshold.solution"];
  out.zerothresholdsearchdir  = control["zero.threshold.searchdir"];
  out.suffdecr                = control["suffdecr.linesearch"];
  out.stepsizereduce          = control["stepsizereduce"];
  out.minstepsize             = control["minstepsize"];
  out.identitycontribincrease = control["identity.contrib.increase"];
  out.maxiteractiveset        = control["maxiter.activeset"];
  out.e                       = control["eps"];
  return out;
}

// Compute a maximum-likelihood estimate (MLE) of the mixture
// proportions in the multinomial mixture model by iterating the
// mix-SQP updates for a fixed number of iterations.
//
// Input argument L is an n x m matrix with non-negative entries;
// input w is a vector of length n containing a non-negative "weight"
// associated with each row of L; input argument x0 is the initial
// estimate of the mixture proportions; input P is a matrix of the
// same dimension as L, and is used to store the posterior mixture
// assignment probabilities; and input "numiter" specifies the number
// of mix-SQP updates to perform.
//
// The return value is a vector of length m containing the updated
// mixture proportions.
//
// Note that L need not be normalized; it will automatically be
// normalized inside this function. Input x0 also need not be
// normalized.
//
// Also note that it does not make sense to compute a MLE of the
// mixture proportions when n < 2 and/or when m < 2; mixsqp will supply
// a result in such cases, but the result will not be valid.
vec mixsqp (const mat& L, const vec& w, const vec& x0, unsigned int numiter,
	    const mixsqp_control_params& control, vec& objective) {
  int m  = L.n_cols;
  mat L1 = L;
  mat Z  = L;
  vec x  = x0;
  mat H(m,m);
  normalizecols(L1);
  mixsqp(L1,w,x,Z,H,numiter,control,objective);
  return x;
}

// Use this variant of mixsqp if you plan on using the same L matrix
// multiple times, or for calling mixsqp multiple times with matrices
// of the same dimension. In the first case, you can reuse the L, Z
// and H matrices; in the latter case, you can reuse the Z and H
// matrices.
//
// For the result to be valid, the matrix L1 should be normalized
// beforehand so that each column sums to 1; Z should be a matrix of
// the same size as L1; and H should be an m x m matrix, where m is the
// number of columns in L1.
//
// Note that input x does not need to be normalized.
void mixsqp (const mat& L1, const vec& w, vec& x, mat& Z, mat& H, 
	     unsigned int numiter, const mixsqp_control_params& control, 
	     vec& objective) {
  double e = control.e;
  int    m = x.n_elem;

  // Initialize storage for matrices and vectors used in the
  // computations below.
  uvec j(m);
  vec  g(m);
  vec  ghat(m);
  vec  y(m);
  vec  xnew(m);
      
  // Normalize the "weights".
  vec w1 = w/sum(w);
  
  // Repeat until the convergence criterion is met, or until we reach
  // the maximum number of (outer loop) iterations.
  for (unsigned int i = 0; i < numiter; i++) {

    // Zero any co-ordinates that are below the specified threshold.
    j = find(x <= control.zerothresholdsolution);
    x(j).fill(0);
    
    // Compute the value of the objective at x.
    objective(i) = compute_objective(L1,w1,x,e);

    // Compute the gradient and Hessian.
    compute_grad(L1,w1,x,e,g,H,Z);
    
    // Solve the quadratic subproblem to obtain a search direction.
    ghat   = g - H*x;
    y      = x;
    activesetqp(H,ghat,y,control.maxiteractiveset,control.
		zerothresholdsearchdir,control.convtolactiveset,
		control.identitycontribincrease);
    
    // Run backtracking line search.
    backtracking_line_search(objective(i),L1,w1,g,x,y,e,control.suffdecr,
			     control.stepsizereduce,control.minstepsize,
			     xnew);
    
    // Update the solution.
    x = xnew;
  }
}

// Return a or b, which ever is smaller.
inline double min (double a, double b) {
  double y;
  if (a < b)
    y = a;
  else
    y = b;
  return y;
}

// Compute the value of the ("modified") objective at x.
double compute_objective (const mat& L, const vec& w, const vec& x, double e) {
  vec u = L*x + e;
  if (u.min() <= 0)
    stop("Objective is -Inf");
  return sum(x) - sum(w % log(u));
}

// Compute the gradient and Hessian of the ("modified") objective at x.
void compute_grad (const mat& L, const vec& w, const vec& x, const double e,
		   vec& g, mat& H, mat& Z) {
  vec u = L*x + e;
  g = -trans(L) * (w/u) + 1;
  Z = L;
  Z.each_col() %= (sqrt(w)/u);
  H = trans(Z) * Z;
}

// Return the largest step size maintaining feasibility (x >= 0) for
// the given the search direction (p);
inline void feasible_stepsize (const vec& x, const vec& p, int& j, double& a) {
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

// This implements the active-set method from p. 472 of of Nocedal &
// Wright, Numerical Optimization, 2nd ed, 2006.
void activesetqp (const mat& H, const vec& g, vec& y, int maxiter,
		  double zerosearchdir, double tol, double ainc) {
  int    m = g.n_elem;
  int    k;
  double a;
  vec    b(m);
  vec    p(m);
  vec    bs(m);
  vec    ps(m);
  mat    Hs(m,m);
  mat    Bs(m,m);
  uvec   i(m);
  uvec   j(m);
  bool   add_to_working_set;

  // This vector is used to keep track of the working set; all zero
  // entries of "t" are co-ordinates belonging to the working set.
  uvec t = (y > 0);
  
  // Run active set method to solve the quadratic subproblem.
  for (int iter = 0; iter < maxiter; iter++) {

    // Find the co-ordinates inside (j) and outside (i) the working
    // set.
    i = find(t != 0);
    j = find(t == 0);

    // Make sure the co-ordinates in the working set are set to zero.
    y(j).fill(0);
    
    // Define the smaller quadratic subproblem.
    Hs    = H(i,i);
    b     = g;
    b(i) += Hs*y(i);
    bs    = b(i);

    // Solve the quadratic subproblem to obtain a search direction.
    p.fill(0);
    compute_activeset_searchdir(Hs,bs,ps,Bs,ainc);
    p(i) = ps;

    // Reset the step size.
    a = 1;
    add_to_working_set = false;
    
    // Check that the search direction is close to zero.
    if (norm(p,"inf") <= zerosearchdir) {

      // Calculate b for all co-ordinates.
      b = g + H*y;
      
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
      if (k >= 0 && a < 1) {

        // Blocking constraint exists; find and add it to the
        // working set (but only if there are two or more non-zero
        // co-ordinates).
	if (i.n_elem > 1)
	  add_to_working_set = true;
      }
      
      // Move to the new iteration along the search direction. The new
      // iterate must have only positive entries.
      y += a*p;
      j  = find(y < 0);
      y(j).fill(0);

      // Make sure all co-ordinates in the working set are set to
      // zero.
      y(j).fill(0);
      if (add_to_working_set) {
	t(k) = 0;
	y(k) = 0;
      }
    }
  }
}

// Get the initial scalar multiplier for the identity matrix based on
// examining the diagonal entries of the Hessian.
inline double init_hessian_correction (const mat& H, double a0) {
  double d = H.diag().min();
  double a;
  if (d > a0)
    a = 0;
  else
    a = a0 - d;
  return a;
}

// This implements Algorithm 3.3, "Cholesky with added multiple of the
// identity", from Nocedal & Wright, 2nd ed, p. 51.
void compute_activeset_searchdir (const mat& H, const vec& y, vec& p,
				  mat& B, double ainc) {
  double a0   = 1e-15;
  double amax = 1e15;
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
    if (a*ainc > amax)
      break;
    else if (chol(R,B))
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
// the search direction is given by p = y - x.
void backtracking_line_search (double f, const mat& L, const vec& w,
			       const vec& g, const vec& x, const vec& y,
			       double e, double suffdecr, double beta,
			       double amin, vec& xnew) {
  int    k;
  double a, afeas;
  double fnew;
  
  // Determine the largest step size maintaining feasibility; if it is
  // larger than the minimum step size, return the minimum step size
  // that maintains feasibility of the solution. Otherwise, continue
  // to backtracking line search.
  vec p = y - x;
  feasible_stepsize(x,p,k,afeas);
  if (afeas <= amin)
    xnew = afeas*y + (1 - afeas)*x;
  else {

    // Set the initial step size.
    a = min(1,afeas);
    
    // Iteratively reduce the step size until either (1) we can't reduce
    // any more (because we have hit the minimum step size constraint),
    // or (2) the new candidate solution satisfies the "sufficient
    // decrease" condition.
    while (true) {
      xnew = a*y + (1 - a)*x;
      fnew = compute_objective(L,w,xnew,e);

      // Check whether the new candidate solution satisfies the
      // sufficient decrease condition, and remains feasible. If so,
      // accept this candidate solution.
      if ((xnew.min() >= 0) &&
	  (fnew <= f + suffdecr*a*dot(y - x,g)))
        break;
    
      // If we cannot decrease the step size further, terminate the
      // backtracking line search, and set the step size to be the
      // minimum step size.
      else if (a*beta < amin) {

        // We need to terminate backtracking line search because we have
        // arrived at the smallest allowable step size.
        a    = amin;
        xnew = a*y + (1 - a)*x;
        if (xnew.min() < 0) {
  	  a    = 0;
	  xnew = x;
        }
        break;
      }
    
      // The new candidate does not satisfy the sufficient decrease
      // condition, so we need to try again with a smaller step size.
      a *= beta;
    }
  }
}
