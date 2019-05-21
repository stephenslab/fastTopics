// This is included to suppress the warnings from solve() when the
// system is singular or close to singular.
#define ARMA_DONT_PRINT_ERRORS

#include <RcppArmadillo.h>

// This is needed to tell R where to find the additional header files.
//
// [[Rcpp::depends(RcppArmadillo)]]
//

using namespace Rcpp;
using namespace arma;

// FUNCTION DECLARATIONS
// ---------------------
void activesetqp(const mat& H, const vec& g, vec& y, uvec& t, int maxiter,
		 double zerothresholdsearchdir, double tol);
void compute_activeset_searchdir (const mat& H, const vec& y, vec& p, mat& B);
void backtracking_line_search (double f, const mat& L, const vec& w,
			       const vec& g, const vec& x, const vec& p,
			       const vec& e, double suffdecr, double beta,
			       double amin, vec& y, vec& u);
void   compute_grad (const mat& L, const vec& w, const vec& x,const vec& e,
		     vec& g, mat& H, mat& Z);
double init_hessian_correction (const mat& H, double a0);
double feasible_stepsize (const vec& x, const vec& p);
double compute_objective (const mat& L, const vec& w, const vec& x,
			  const vec& e, vec& u);
double min (double a, double b);

// FUNCTION DEFINITIONS
// --------------------
// SQP algorithm for computing a maximum-likelihood estimate of a
// mixture model. For more information, see the help and comments
// accompanying the mixsqp R function and, in particular, see how
// mixsqp_rcpp is called inside the mixsqp function.
// 
// [[Rcpp::export]]
vec mixsqp_rcpp (const mat& L, const vec& w, const vec& x0, 
		 double tol, double zerothreshold,
		 double zerosearchdir, double suffdecr,
		 double stepsizereduce, double minstepsize,
		 const vec& e, int numiter, int maxiteractiveset,
		 bool verbose) {
  
  // Get the number of rows (n) and columns (m) of the conditional
  // likelihood matrix.
  int n = L.n_rows;
  int m = L.n_cols;

  // Initialize the solution.
  vec x = x0;
  
  // Scalars, vectors and matrices used in the computations below.
  double obj;
  uvec i(m);
  vec  g(m);
  vec  ghat(m);
  vec  p(m);
  vec  u(n);
  mat  H(m,m);
  mat  Z(n,m);
  uvec t(m);
  vec  y(m);
  
  // Iterate the SQP updates for a fixed number of iterations.
  for (int iter = 0; iter < numiter; iter++) {

    // Zero any co-ordinates that are below the specified threshold.
    i = find(x <= zerothreshold);
    x.elem(i).fill(0);
    
    // Compute the value of the objective at x.
    obj = compute_objective(L,w,x,e,u);
    if (verbose)
      Rprintf("%4d %+0.15f\n",iter,obj);

    // Compute the gradient and Hessian.
    compute_grad(L,w,x,e,g,H,Z);
    
    // This specifies the "inactive set".
    t = (x > 0);

    // Solve the quadratic subproblem to obtain a search direction.
    ghat = g - H*x;
    activesetqp(H,ghat,y,t,maxiteractiveset,zerosearchdir,tol);
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
void activesetqp (const mat& H, const vec& g, vec& y, uvec& t,
		  int maxiter, double zerosearchdir, double tol) {
  int    m     = g.n_elem;
  double nnz   = sum(t);
  double alpha;  // The step size.
  int    k;
  vec    b(m);   // Vector of length m storing H*y + 2*g + 1.
  vec    p(m);   // Vector of length m storing the search direction.
  vec    p0(m);  // Search direction for nonzero co-ordinates only.
  vec    bs(m);
  vec    z(m);
  mat    Hs(m,m);
  mat    B(m,m);
  uvec   S(m);
  uvec   i0(m);
  uvec   i1(m);
  
  // Initialize the solution to the QP subproblem, y.
  y.fill(0);
  i1 = find(t);
  y.elem(i1).fill(1/nnz);
    
  // Run active set method to solve the QP subproblem.
  for (int j = 0; j < maxiter; j++) {
    
    // Define the smaller QP subproblem.
    i0 = find(1 - t);
    i1 = find(t);
    b  = H*y + g;
    bs = b.elem(i1);
    Hs = H.elem(i1,i1);
      
    // Solve the smaller problem.
    p.fill(0);
    compute_activeset_searchdir(Hs,bs,p0,B);
    p.elem(i1) = p0;
      
    // Reset the step size.
    alpha = 1;
    
    // Check that the search direction is close to zero.
    if ((p.max() <= zerosearchdir) & (-p.min() <= zerosearchdir)) {
        
      // If all the Lagrange multiplers in the working set (that is,
      // zeroed co-ordinates) are positive, or nearly positive, we
      // have reached a suitable solution.
      if (i0.n_elem == 0) {
	j++;
	break;
      } else if (b(i0).min() >= -tol) {
	j++;
	break;
      }

      // Find an co-ordinate with the smallest multiplier, and remove
      // it from the working set.
      k    = i0[b(i0).index_min()];
      t[k] = 1;

    // In this next part, we consider adding a co-ordinate to the
    // working set (but only if there are two or more non-zero
    // co-ordinates).
    } else {
        
      // Define the step size.
      alpha = 1;
      p0    = p;
      p0.elem(i0).fill(0);
      S = find(p0 < 0);
      if (!S.is_empty()) {
        z = -y.elem(S)/p.elem(S);
        k = z.index_min();
        if (z[k] < 1) {
            
          // Blocking constraint exists; find and add it to the
          // working set (but only if there are two or more non-zero
          // co-ordinates).
          alpha = z[k];
	  if (i1.n_elem >= 2)
	    t[S[k]] = 0;
        }
      }
      
      // Move to the new "inner loop" iterate (y) along the search
      // direction.
      y += alpha * p;
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
  double fnew;

  // Determine the largest step size maintaining feasibility; if it is
  // larger than the minimum step size, return the minimum step size
  // that maintains feasibility of the solution. Otherwise, continue
  // to backtracking line search.
  double afeas = feasible_stepsize(x,p);
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
double feasible_stepsize (const vec& x, const vec& p) {
  uvec i = find(p < 0);
  vec  t = -x.elem(i)/p.elem(i);
  return t.min();
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

// Return a or b, which ever is smaller.
double min (double a, double b) {
  double y;
  if (a < b)
    y = a;
  else
    y = b;
  return y;
}
