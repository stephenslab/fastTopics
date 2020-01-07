#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// FUNCTION DECLARATIONS
// ---------------------
double compute_objective (const mat& L, const vec& w, const vec& x,
			  const vec& e);
void   compute_grad      (const mat& L, const vec& w, const vec& x,
			  const vec& e, vec& g, mat& H, mat& Z);
void   activesetqp       (const mat& H, const vec& g, vec& y, int maxiter,
			  double zerosearchdir, double tol, double ainc);
void   compute_activeset_searchdir (const mat& H, const vec& y, vec& p, mat& B,
				    double ainc);
void   backtracking_line_search (double f, const mat& L, const vec& w,
				 const vec& g, const vec& x, const vec& y,
				 const vec& e, double suffdecr, double beta,
				 double amin, vec& xnew);

// FUNCTION DEFINITIONS
// --------------------
// SQP algorithm for computing a maximum-likelihood estimate of a
// mixture model. For more information, see the help and comments
// accompanying the mixsqp R function and, in particular, see how
// mixsqp_rcpp is called inside the mixsqp function.
// 
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List mixsqp_rcpp (const arma::mat& L, const arma::vec& w, const arma::vec& x0,
		  double convtolactiveset, double zerothresholdsolution,
		  double zerothresholdsearchdir, double suffdecr,
		  double stepsizereduce, double minstepsize,
		  double identitycontribincrease, const arma::vec& eps,
		  int numitersqp, int maxiteractiveset) {
  
  // Get the number columns of the data matrix.
  int m = L.n_cols;

  // Initialize storage for the "obj" output.
  vec obj(numitersqp,fill::zeros);

  // Initialize the solution.
  vec x = x0;
  
  // Initialize storage for matrices and vectors used in the
  // computations below.
  vec  g(m);
  vec  ghat(m);
  mat  H(m,m);
  uvec j(m);
  vec  y(m);
  vec  d(m);
  vec  xnew(m);
  mat  Z = L;
      
  // Repeat until the convergence criterion is met, or until we reach
  // the maximum number of (outer loop) iterations.
  for (uint i = 0; i < numitersqp; i++) {

    // Zero any co-ordinates that are below the specified threshold.
    j = find(x <= zerothresholdsolution);
    x(j).fill(0);
    
    // Compute the value of the objective at x.
    obj(i) = compute_objective(L,w,x,eps);

    // Compute the gradient and Hessian.
    compute_grad(L,w,x,eps,g,H,Z);
    
    // Solve the quadratic subproblem to obtain a search direction.
    ghat   = g - H*x + 1;
    y      = x;
    activesetqp(H,ghat,y,maxiteractiveset,zerothresholdsearchdir,
		convtolactiveset,identitycontribincrease);
    
    // Run backtracking line search.
    backtracking_line_search(obj(i),L,w,g,x,y,eps,suffdecr,stepsizereduce,
			     minstepsize,xnew);
    
    // Update the solution.
    x = xnew;
  }

  // CONSTRUCT OUTPUT
  // ----------------
  return List::create(Named("x") = x,Named("objective") = obj);
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

// Compute the value of the (unmodified) objective at x.
double compute_objective (const mat& L, const vec& w, const vec& x,
			  const vec& e) {
  vec u = L*x + e;
  if (u.min() <= 0)
    stop("Objective is -Inf");
  return -sum(w % log(u));
}

// Compute the gradient and Hessian of the (unmodified) objective at x.
void compute_grad (const mat& L, const vec& w, const vec& x, const vec& e,
		   vec& g, mat& H, mat& Z) {
  vec u = L*x + e;
  g = -trans(L) * (w/u);
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
    if ((p.max() <= zerosearchdir) && (-p.min() <= zerosearchdir)) {

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
  double amax = 1;
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
    if (chol(R,B) || (a*ainc > amax))
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
// the search direction is given by p = y - x.
void backtracking_line_search (double f, const mat& L, const vec& w,
			       const vec& g, const vec& x, const vec& y,
			       const vec& e, double suffdecr, double beta,
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
  	 (fnew + sum(xnew) <= f + sum(x) + suffdecr*a*dot(y - x,g + 1)))
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
