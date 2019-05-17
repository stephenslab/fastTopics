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

// FUNCTION DECLARATIONS
// ---------------------
double mixobjective (const mat& L, const vec& w, const vec& x, const vec& e, 
		     vec& u);
void   computegrad  (const mat& L, const vec& w, const vec& x, const vec& e, 
		     vec& g, mat& H, vec& u, mat& Z);
double activesetqp  (const mat& H, const vec& g, vec& y, uvec& t,
		     int maxiteractiveset, double zerothresholdsearchdir, 
		     double convtolactiveset, double identitycontribincrease);
void computeactivesetsearchdir (const mat& H, const vec& y, vec& p,
				mat& B, double ainc);
void backtrackinglinesearch (double f, const mat& L, const vec& w,
			     const vec& g, const vec& x, const vec& p,
			     const vec& eps, double suffdecr,
			     double stepsizereduce, double minstepsize, 
			     double& nls, double& stepsize, vec& y, vec& u);

// FUNCTION DEFINITIONS
// --------------------
// SQP algorithm for computing a maximum-likelihood estimate of a
// mixture model. For more information, see the help and comments
// accompanying the mixsqp R function and, in particular, see how
// mixsqp_rcpp is called inside the mixsqp function.
// 
// [[Rcpp::export]]
List mixsqp_rcpp (const arma::mat& L, const arma::vec& w, const arma::vec& x0, 
                  double convtolsqp, double convtolactiveset,
		  double zerothresholdsolution, double zerothresholdsearchdir,
		  double suffdecr, double stepsizereduce, double minstepsize,
		  double identitycontribincrease, const arma::vec& eps,
		  int maxitersqp, int maxiteractiveset, bool verbose) {
  
  // Get the number of rows (n) and columns (m) of the conditional
  // likelihood matrix.
  int n = L.n_rows;
  int m = L.n_cols;

  // PREPARE DATA STRUCTURES
  // -----------------------
  // Initialize storage for the outputs obj, gmin, nnz, nqp and dmax.
  vec obj(maxitersqp);
  vec gmin(maxitersqp);
  vec nnz(maxitersqp);
  vec nqp(maxitersqp);
  vec nls(maxitersqp);
  vec stepsize(maxitersqp);
  vec dmax(maxitersqp);
  obj.zeros();
  gmin.zeros();
  nnz.zeros();
  nqp.fill(-1);
  nls.fill(-1);
  stepsize.fill(-1);
  dmax.fill(-1);

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
  
  double status = 1;  // Convergence status.
  
  // Initialize some loop variables used in the loops below.
  int i = 0; 
  
  // Repeat until the convergence criterion is met, or until we reach
  // the maximum number of (outer loop) iterations.
  for (i = 0; i < maxitersqp; i++) {

    // Compute the value of the objective at x (obj).
    obj[i] = mixobjective(L,w,x,eps,u);

    // COMPUTE GRADIENT AND HESSIAN
    // ----------------------------
    computegrad(L,w,x,eps,g,H,u,Z);
    
    // Determine the nonzero co-ordinates in the current estimate of
    // the solution, x. This specifies the "inactive set".
    t = (x >= zerothresholdsolution);

    // Report on the algorithm's progress. Here we compute: the
    // smallest gradient value (gmin), which is used as a convergence
    // criterion; and the number of nonzeros in the solution (nnz).
    // Note that only the dual residuals (gmin's) corresponding to the
    // nonzero co-ordinates are relevant.
    gmin[i] = 1 + g.min();
    nnz[i]  = sum(t);
    if (verbose) {
      if (i == 0)
        Rprintf("%4d %+0.9e %+0.3e%4d  ------   ------   --  --\n",
		i + 1,obj[i],-gmin[i],int(nnz[i]));
      else
        Rprintf("%4d %+0.9e %+0.3e%4d %0.2e %0.2e %3d %3d\n",i + 1,obj[i],
		-gmin[i],(int) nnz[i],stepsize[i-1],dmax[i-1],
		(int) nqp[i-1],(int) nls[i-1]);
    }
    
    // CHECK CONVERGENCE
    // -----------------
    // Convergence is reached with the maximum dual residual is
    // small. The negative of "gmin" is also the maximum dual residual
    // (denoted as "rdual" on p. 609 of Boyd & Vandenberghe, "Convex
    // Optimization", 2009). Although "gmin" here includes both zero
    // (active) and non-zero (inactive) co-ordinates, this condition
    // is trivially satisfied for the zero co-ordinates as the
    // gradient must be non-negative for these co-ordinates. (See
    // communication with Youngseok on Slack.)
    if (gmin[i] >= -convtolsqp) {
      status = 0;
      i++;
      break;
    }

    // This is also a good point to check for a user interrupt; if the
    // user requests an interrupt, then an exception is thrown and
    // control is returned to the R console.
    Rcpp::checkUserInterrupt();
    
    // SOLVE QUADRATIC SUBPROBLEM
    // --------------------------
    // Run the active-set solver to obtain a search direction.
    ghat   = g - H*x + 1;
    nqp[i] = activesetqp(H,ghat,y,t,maxiteractiveset,zerothresholdsearchdir,
			 convtolactiveset,identitycontribincrease);
    p = y - x;
    
    // BACKTRACKING LINE SEARCH
    // ------------------------
    backtrackinglinesearch(obj[i],L,w,g,x,p,eps,suffdecr,stepsizereduce,
			   minstepsize,nls[i],stepsize[i],y,u);
    
    // UPDATE THE SOLUTION
    // -------------------
    d       = abs(x - y);
    dmax[i] = d.max();
    x = y;
  }

  // CONSTRUCT OUTPUT
  // ----------------
  return List::create(Named("x")         = x,
		      Named("status")    = status,
		      Named("objective") = obj.head(i),
		      Named("max.rdual") = -gmin.head(i),
		      Named("nnz")       = nnz.head(i),
		      Named("stepsize")  = stepsize.head(i),
		      Named("max.diff")  = dmax.head(i),
		      Named("nqp")       = nqp.head(i),
		      Named("nls")       = nls.head(i));
}

// Compute the value of the (unmodified) objective at x, assuming x is
// (primal) feasible; arguments L and w specify the objective, and e
// is an additional positive constant near zero. Input argument u is a
// vector of length n used to store an intermediate result used in the
// calculation of the objective.
double mixobjective (const mat& L, const vec& w, const vec& x, const vec& e,
		     vec& u) {
  u = L*x + e;
  if (u.min() <= 0)
    Rcpp::stop("Halting because the objective function has a non-finite value (logarithms of numbers less than or equal to zero) at the current estimate of the solution");
  return -sum(w % log(u));
}

// Compute the gradient and Hessian of the (unmodified) objective at
// x, assuming x is (primal) feasible; arguments L and w specify the
// objective, and e is an additional positive constant near
// zero. Inputs g and H store the value of the gradient and Hessian (a
// vector of length m, and an m x m matrix). Intermediate results used
// in these calculations are stored in three variables: u, a vector of
// length n; Z, an n x m matrix; and ZW, another n x m matrix.
void computegrad (const mat& L, const vec& w, const vec& x, const vec& e,
		  vec& g, mat& H, vec& u, mat& Z) {
   
  // Compute the gradient g = -L'*(w.*u) where u = 1./(L*x + e), ".*"
  // denotes element-wise multiplication, and "./" denotes
  // element-wise division.
  u = L*x + e;
  u = w / u;
  g = -trans(L) * u;
  
  // Compute the Hessian H = L'*U*W*U*L, where U = diag(u) and 
  // W = diag(w), with vector u defined as above.
  Z = L;
  u /= sqrt(w);
  Z.each_col() %= u;
  H = trans(Z) * Z;
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
			     double& nls, double& stepsize, vec& y, vec& u) {
  double fnew;
  double newstepsize;
  stepsize = 0.99;
  nls      = 0;

  // Iteratively reduce the step size until either (1) we can't reduce
  // any more (because we have hit the minimum step size constraint),
  // or (2) the new candidate solution satisfies the "sufficient
  // decrease" condition.
  while (true) {
    y    = x + stepsize*p;
    fnew = mixobjective(L,w,y,eps,u);
    nls++;

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
