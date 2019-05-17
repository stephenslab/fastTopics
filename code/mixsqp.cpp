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
double activesetqp  (const mat& H, const vec& g, vec& y, uvec& t,
		     int maxiteractiveset, double zerothresholdsearchdir, 
		     double convtolactiveset, double identitycontribincrease);
void computeactivesetsearchdir (const mat& H, const vec& y, vec& p,
				mat& B, double ainc);
void backtrackinglinesearch (double f, const mat& L, const vec& w,
			     const vec& g, const vec& x, const vec& p,
			     const vec& eps, double suffdecr,
			     double stepsizereduce, double minstepsize, 
			     vec& y, vec& u);
double mixobjective (const mat& L, const vec& w, const vec& x,
		     const vec& e, vec& u);
void   computegrad  (const mat& L, const vec& w, const vec& x,
		     const vec& e, vec& g, mat& H, mat& Z);

// FUNCTION DEFINITIONS
// --------------------
// SQP algorithm for computing a maximum-likelihood estimate of a
// mixture model. For more information, see the help and comments
// accompanying the mixsqp R function and, in particular, see how
// mixsqp_rcpp is called inside the mixsqp function.
// 
// [[Rcpp::export]]
vec mixsqp_rcpp (const mat& L, const vec& w, const vec& x0, 
		 double convtolactiveset,
		 double zerothresholdsolution, double zerothresholdsearchdir,
		 double suffdecr, double stepsizereduce, double minstepsize,
		 double identitycontribincrease, const vec& eps,
		 int numiter, int maxiteractiveset) {
  
  // Get the number of rows (n) and columns (m) of the conditional
  // likelihood matrix.
  int n = L.n_rows;
  int m = L.n_cols;

  // Initialize the solution.
  vec x = x0;
  
  // Initialize storage for matrices and vectors used in the
  // computations below.
  vec  g(m);    // Vector of length m storing the gradient.
  vec  ghat(m); // Vector of length m storing gradient of subproblem.
  vec  p(m);    // Vector of length m storing the search direction.
  vec  u(n);    // Vector of length n storing L*x + eps or its log.
  mat  H(m,m);  // m x m matrix storing Hessian.
  mat  Z(n,m);  // n x m matrix Z = D*L, where D = diag(1/(L*x+e)).
  uvec t(m);    // Temporary unsigned int. vector result of length m.
  vec  y(m);    // Vector of length m storing y
  vec  d(m);    // Vector of length m storing absolute
                // differences between between two solution
	        // estimates.
  
  // Iterate the SQP updates for a fixed number of iterations.
  for (int iter = 0; iter < numiter; iter++) {

    // Compute the value of the objective at x.
    double obj = mixobjective(L,w,x,eps,u);

    // *** FOR TESTING ***
    Rprintf("%+0.12f\n",obj);

    // Compute the gradient and Hessian.
    computegrad(L,w,x,eps,g,H,Z);
    
    // Determine the nonzero co-ordinates in the current estimate of
    // the solution, x. This specifies the "inactive set".
    t = (x >= zerothresholdsolution);

    // Solve the quadratic subproblem to obtain a search direction.
    ghat = g - H*x;
    activesetqp(H,ghat,y,t,maxiteractiveset,zerothresholdsearchdir,
		convtolactiveset,identitycontribincrease);
    p = y - x;
    
    // Run backtracking line search.
    backtrackinglinesearch(obj,L,w,g,x,p,eps,suffdecr,stepsizereduce,
			   minstepsize,y,u);
    x = y;
  }

  return x;
}

// This implements the active-set method from p. 472 of of Nocedal &
// Wright, Numerical Optimization, 2nd ed, 2006.
double activesetqp (const mat& H, const vec& g, vec& y, uvec& t,
		    int maxiteractiveset, double zerothresholdsearchdir, 
		    double convtolactiveset, double identitycontribincrease) {
  int    m     = g.n_elem;
  double nnz   = sum(t);
  double alpha;  // The step size.
  int    j, k;
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
  for (j = 0; j < maxiteractiveset; j++) {
    
    // Define the smaller QP subproblem.
    i0 = find(1 - t);
    i1 = find(t);
    b  = H*y + g;
    bs = b.elem(i1);
    Hs = H.elem(i1,i1);
      
    // Solve the smaller problem.
    p.fill(0);
    computeactivesetsearchdir(Hs,bs,p0,B,identitycontribincrease);
    p.elem(i1) = p0;
      
    // Reset the step size.
    alpha = 1;
    
    // Check that the search direction is close to zero (according to
    // the "zerothresholdsearchdir" parameter).
    if ((p.max() <= zerothresholdsearchdir) &
	(-p.min() <= zerothresholdsearchdir)) {
        
      // If all the Lagrange multiplers in the working set (that is,
      // zeroed co-ordinates) are positive, or nearly positive, we
      // have reached a suitable solution.
      if (i0.n_elem == 0) {
	j++;
	break;
      } else if (b(i0).min() >= -convtolactiveset) {
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

  return j;
}

// This implements Algorithm 3.3, "Cholesky with added multiple of the
// identity", from Nocedal & Wright, 2nd ed, p. 51.
void computeactivesetsearchdir (const mat& H, const vec& y, vec& p,
				mat& B, double ainc) {
  double a0   = 1e-15;
  double amax = 1;
  int    n = y.n_elem;
  double d = H.diag().min();
  double a = 0;
  mat    I(n,n,fill::eye);
  mat    R(n,n);

  // Initialize the scalar multiplier for the identity matrix.
  if (d < 0)
    a = a0 - d;

  // Repeat until a modified Hessian is found that is symmetric
  // positive definite, or until we cannot modify it any longer.
  while (true) {

    // Compute the modified Hessian.
    B = H + a*I;
    
    // Attempt to compute the Cholesky factorization of the modified
    // Hessian. If this fails, increase the contribution of the
    // identity matrix in the modified Hessian.
    if (chol(R,B) | (a*ainc > amax))
      break;
    else if (a <= 0)
      a = a0;
    else
      a *= ainc;
  }
  
  // Compute the search direction using the modified Hessian.
  p = -solve(B,y);
}

// This implements the backtracking line search algorithm from p. 37
// of Nocedal & Wright, Numerical Optimization, 2nd ed, 2006. Note
// that sum(x) = sum(y) = 1, so replacing g by g + 1 in dot product of
// p and g has no effect.
void backtrackinglinesearch (double f, const mat& L, const vec& w,
			     const vec& g, const vec& x, const vec& p,
			     const vec& eps, double suffdecr,
			     double stepsizereduce, double minstepsize, 
			     vec& y, vec& u) {
  double fnew;
  double newstepsize;
  double stepsize = 0.99;

  // Iteratively reduce the step size until either (1) we can't reduce
  // any more (because we have hit the minimum step size constraint),
  // or (2) the new candidate solution satisfies the "sufficient
  // decrease" condition.
  while (true) {
    y    = x + stepsize*p;
    fnew = mixobjective(L,w,y,eps,u);

    // Check whether the new candidate solution (y) satisfies the
    // sufficient decrease condition, and remains feasible. If so,
    // accept this candidate solution.
    if (y.min() >= 0 & fnew <= f + suffdecr*stepsize*dot(p,g))
      break;
    newstepsize = stepsizereduce * stepsize;
    if (newstepsize < minstepsize) {

      // We need to terminate backtracking line search because we have
      // arrived at the smallest allowed step size.
      stepsize = minstepsize;
      y        = x + stepsize*p;
      if (y.min() < 0) {
	stepsize = 0;
	y        = x;
      }
      break;
    }
    
    // The new candidate does not satisfy the sufficient decrease
    // condition, so we need to try again with a smaller step size.
    stepsize = newstepsize;
  }
}

// Compute the value of the objective at x; arguments L and w specify
// the objective, and e is a vector in which the entries can be set to
// small, positive numbers, or to zero. Input u stores an intermediate
// result used in the calculation.
double mixobjective (const mat& L, const vec& w, const vec& x, const vec& e, 
		     vec& u) {
  u = L*x + e;
  if (u.min() <= 0)
    Rcpp::stop("Objective is -Inf");
  return sum(x) - sum(w % log(u));
}

// Compute the gradient and Hessian of the objective at x.
void computegrad (const mat& L, const vec& w, const vec& x,
		  const vec& e, vec& g, mat& H, mat& Z) {
  vec u = L*x + e;
  g = -trans(L)*(w/u) + 1;
  Z = L;
  Z.each_col() %= (sqrt(w)/u);
  H = trans(Z) * Z;
}
