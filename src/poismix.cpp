#include "poismix.h"
#include "misc.h"
#include "mixem.h"

using namespace arma;

// INLINE FUNCTION DEFINITIONS
// ---------------------------
// This is the same as compute_grad_scd (see below), except that u =
// sum(L[,j]), w must contain only the nonzero values, and vectors l
// and Lx should contain only the entries associated with these
// nonzero values (vectors l, w and Lx should be the same length). 
// This is used by scd_kl_update.
inline void compute_grad_scd_sparse (const vec& l, const vec& w, const vec& Lx,
				     double u, double x, double& g, double& h,
				     vec& r, double e) {
  r = l/(Lx + e);
  h = dot(square(r),w);
  g = dot(r,w) - u + h*x;
}

// Compute the first-order and second-order partial derivatives of the
// loss function (i.e., the negative Poisson log-likelihood) with
// respect to the jth mixture weight, x[j], where: input argument l =
// L[,j] is the jth column of the n x m data matrix, L; input w is a
// vector of length n containing a count or "pseudocount" associated
// with each row of L; x is the initial estimate of the mixture
// weights; and Lx = L*x. Input argument u is a vector storing an
// intermediate result of the same size as w. The outputs are g = h*x
// - d and h, where d and h are the first-order and second-order
// partial derivatives of the loss function. This is used by
// scd_kl_update.
inline void compute_grad_scd (const vec& l, const vec& w, const vec& Lx,
			      double x, double& g, double& h, vec& r,
			      double e) {
  compute_grad_scd_sparse(l,w,Lx,sum(l),x,g,h,r,e);
}

// Given g = h*x - d and h, where d and h are the first-order and
// second-order partial derivatives of the loss function (negative
// Poisson log-likelihood), compute the SCD update of the mixture
// weight, x[j]. This is used by scd_kl_update.
inline double project_iterate_scd (double g, double h, double e) {
  double x = g/(h + e);
  x = maximum(x,0);
  return x;
}

// Compute the first-order (g) and second-order (h) partial
// derivatives of the loss function (i.e., the negative Poisson
// log-likelihood) with respect to the jth mixture weight, x[j]. See
// compute_grad_scd for an explanation of the inputs. This is used by
// ccd_kl_update.
inline void compute_grad_ccd (const vec& l, const vec& w, const vec& Lx,
			      double x, double& g, double& h, double e) {
  unsigned int n = l.n_elem;
  double       t;
  g = 0;
  h = 0;
  for (unsigned int i = 0; i < n; i++) {
    t  = w(i)/(Lx(i) + e);
    g += l(i)*(1 - t);
    h += l(i)*l(i)*t/(Lx(i) + e);
  }
}

// This is the same as compute_grad_ccd, except that u = sum(L[,j]), w
// must contain only the nonzero values, and vectors l and Lx should
// contain only the entries associated with these nonzero values
// (vectors l, w and Lx should be the same length). This is used by
// ccd_kl_update.
inline void compute_grad_ccd_sparse (const vec& l, const vec& w, const vec& Lx,
				     double u, double x, double& g, double& h,
				     double e) {
  unsigned int n = w.n_elem;
  double       t;
  g = u;
  h = 0;
  for (unsigned int i = 0; i < n; i++) {
    t  = w(i)/(Lx(i) + e);
    g -= l(i)*t;
    h += l(i)*l(i)*t/(Lx(i) + e);
  }
}

// Given the first-order (g) and second-order (h) partial derivatives
// of the loss function (negative Poisson log-likelihood), compute the
// CCD update of the mixture weight, w[j]. This is used by ccd_kl_update.
inline double project_iterate_ccd (double x, double g, double h, double e) {
  double y;
  y = x - g/h + e;
  y = maximum(y,e);
  return y;
}

// FUNCTION DEFINITIONS
// --------------------
// This is mainly used to test the first variant of the poismixem C++
// function.
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec poismixem_rcpp (const arma::mat& L, const arma::vec& w,
			  const arma::vec& x0, unsigned int numiter) {
  return poismixem(L,w,x0,numiter);
}

// This is mainly used to test the second variant of the poismixem C++
// interface.
//
// [[Rcpp::export]]
arma::vec poismixem2_rcpp (const arma::mat& L1, const arma::vec& w,
			   const arma::vec& u, const arma::vec& x0,
			   unsigned int numiter) {
  vec x = x0;
  mat P = L1;
  poismixem(L1,u,w,x,P,numiter);
  return x;
}

// This is mainly used to test the third variant of the poismixem C++
// interface.
//
// [[Rcpp::export]]
arma::vec poismixem3_rcpp (const arma::mat& L1, const arma::vec& w,
			   const arma::vec& u, const arma::uvec& i,
			   const arma::vec& x0, unsigned int numiter) {
  vec x = x0;
  poismixem(L1,u,w,i,x,numiter);
  return x;
}

// This is mainly used to test the first variant of the C++ function
// scd_kl_update.
//
// [[Rcpp::export]]
arma::vec scd_kl_update_rcpp (const arma::mat& L, const arma::vec& w, 
			      const arma::vec& x0, unsigned int numiter, double e) {
  return scd_kl_update(L,w,x0,numiter,e);
}

// This is mainly used to test the second variant of the C++ function
// scd_kl_update.
//
// [[Rcpp::export]]
arma::vec scd_kl_update2_rcpp (const arma::mat& L, const arma::vec& u,
			       const arma::vec& w, const arma::vec& x0,
			       unsigned int numiter, double e) {
  return scd_kl_update(L,u,w,x0,numiter,e);
}

// This is mainly used to test the first variant of the C++ function
// ccd_kl_update.
//
// [[Rcpp::export]]
arma::vec ccd_kl_update_rcpp (const arma::mat& L, const arma::vec& w,
			      const arma::vec& x0, unsigned int numiter, 
			      double e) {
  vec x = x0;
  for (unsigned int iter = 0; iter < numiter; iter++)
    x = ccd_kl_update(L,w,x,e);
  return x;
}

// This is mainly used to test the second variant of the C++ function
// ccd_kl_update.
//
// [[Rcpp::export]]
arma::vec ccd_kl_update2_rcpp (const arma::mat& L, const arma::vec& u,
			       const arma::vec& w, const arma::vec& x0,
			       unsigned int numiter, double e) {
  vec x = x0;
  for (unsigned int iter = 0; iter < numiter; iter++) 
    x = ccd_kl_update(L,u,w,x,e);
  return x;
}

// Compute a maximum-likelihood estimate (MLE) of the mixture weights
// in a Poisson mixture model by iterating the multinomial mixture
// model EM updates for a fixed number of iterations.
//
// Input argument L is an n x m matrix with non-negative entries;
// input w is a vector of length n containing a count or "pseudocount"
// associated with each row of L; input argument x0 is the initial
// estimate of the mixture weights; and input "numiter" specifies the
// number of EM updates to perform.
//
// The return value is a vector of length m containing the updated
// mixture weights.
vec poismixem (const mat& L, const vec& w, const vec& x0, 
	       unsigned int numiter) {
  mat L1 = L;
  mat P  = L;
  vec u  = trans(sum(L,0));
  vec x  = x0;
  normalizecols(L1);
  poismixem(L1,u,w,x,P,numiter);
  return x;
}

// Use this variant of poismixem if you plan on calling poismixem
// multiple times with the same matrix, L. In this call, input u should
// contain the column sums, u = colSums(L), in which L is the matrix
// prior to normalization, and matrix L1 is the normalized version of
// L in which each column sums to 1; that is, L1 = normalize.cols(L).
// Input P should be a matrix of the same size as L1.
//
// Note that in this variant of poismixem, L1 and w do not need to
// contain all the data; once u is correctly calculated, any rows of
// L1 associated with zero weights have no effect, so only the vector
// of nonzero weights w, and the rows of L1 associated with those
// weights, need to be supplied.
void poismixem (const mat& L1, const vec& u, const vec& w, vec& x, mat& P, 
		unsigned int numiter) {
  double s = sum(w);
  
  // Recover the mixture proportions of the multinomial mixture model
  // from the mixture weights of the Poisson mixture model. 
  x %= u;
  x /= sum(x);

  // Perform one or more EM updates for the multinomial mixture model.
  mixem(L1,w,x,P,numiter);
  
  // Recover the mixture weights of the Poisson mixture model from the
  // mixture weights of the multinomial mixture model.
  x *= s;
  x /= u;
}

// This third poismixem interface is similar to the second, except
// that the indices of the nonzero counts are supplied in vector i,
// and w must contain only the nonzero counts (vectors i and w should
// be the same length).
//
// To illustrate, in this example x1 and x2 should be the same after
// running this R code:
//
//   i  <- which(w > 0)
//   x1 <- poismixem2_rcpp(L1,w,u,x0,numiter)
//   x2 <- poismixem3_rcpp(L1,w[i],u,i-1,x0,numiter)
//
void poismixem (const mat& L1, const vec& u, const vec& w, const uvec& i,
		vec& x, unsigned int numiter) {
  unsigned int n = i.n_elem;
  unsigned int m = x.n_elem;
  mat  P(n,m);
  poismixem(L1.rows(i),u,w,x,P,numiter);
}

// Compute a maximum-likelihood estimate (MLE) of the mixture weights
// in a Poisson mixture model by iterating coordinate-wise updates for
// a fixed number of iterations.
//
// Input argument L is an n x m matrix with non-negative entries;
// input w is a vector of length n containing a count or "pseudocount"
// associated with each row of L; input argument x0 is the initial
// estimate of the mixture weights; and input "numiter" specifies the
// number of updates (full passes) to perform.
//
// The return value is a vector of length m containing the updated
// mixture weights.
//  
// This function implements the core part of the sequential
// co-ordinate descent (SCD) algorithm. This C++ code is adapted from
// the code by Xihui Lin and Paul Boutros, which is available for
// download at https://github.com/linxihui/NNLM.
vec scd_kl_update (const mat& L, const vec& w, const vec& x0,
		   unsigned int numiter, double e) {
  unsigned int n = L.n_rows;
  unsigned int m = L.n_cols;
  vec    x  = x0;
  vec    Lx = L * x;
  vec    l(n);
  vec    r(n);
  double h, g, xj, xjnew;
  for (unsigned int iter = 0; iter < numiter; iter++)
    for (unsigned int j = 0; j < m; j++) {
      l     = L.col(j);
      xj    = x(j);
      compute_grad_scd(l,w,Lx,xj,g,h,r,e);
      xjnew = project_iterate_scd(g,h,e);
      Lx   += (xjnew - xj) * l;
      x(j)  = xjnew;
    }
  return x;
}

// In this alternative implementation of the SCD updates, w must
// contain only the nonzero counts, and L must contain only the rows
// of the matrix associated with the nonzero counts (vectors i and w
// should be the same length, and the number of rows in L should be
// the same as the lengths of vectors i and w).
vec scd_kl_update (const mat& L, const vec& u, const vec& w,
		   const vec& x0, unsigned int numiter, double e) {
  unsigned int n = w.n_elem;
  unsigned int m = L.n_cols;
  vec    x  = x0;
  vec    Lx = L * x;
  vec    l(n);
  vec    r(n);
  double h, g, xj, xjnew;
  for (unsigned int iter = 0; iter < numiter; iter++)
    for (unsigned int j = 0; j < m; j++) {
      l     = L.col(j);
      xj    = x(j);
      compute_grad_scd_sparse(l,w,Lx,u(j),xj,g,h,r,e);
      xjnew = project_iterate_scd(g,h,e);
      Lx   += (xjnew - xj) * l;
      x(j)  = xjnew;
    }
  return x;
}

// Implements a single iteration of the cyclic co-ordinate descent (CCD)
// algorithm. See scd_kl_update for an explanation of the inputs.
//
// This implementation is adapted from the C++ code developed by
// Cho-Jui Hsieh and Inderjit Dhillon, which is available for download
// at www.cs.utexas.edu/~cjhsieh/nmf.
vec ccd_kl_update (const mat& L, const vec& w, const vec& x0, double e) {
  unsigned int n = L.n_rows;
  unsigned int m = L.n_cols;
  vec    x  = x0;
  vec    Lx = L * x;
  vec    l(n);
  double g, h, xj, xjnew;
  for (unsigned int j = 0; j < m; j++) {
    l     = L.col(j);
    xj    = x(j);
    compute_grad_ccd(l,w,Lx,xj,g,h,e);
    xjnew = project_iterate_ccd(xj,g,h,e);
    Lx   += (xjnew - xj) * l;
    x(j)  = xjnew;
  }

  return x;
}

// In this alternative implementation of the CCD updates, w must
// contain only the nonzero counts, and L must contain only the rows
// of the matrix associated with the nonzero counts (vectors i and w
// should be the same length, and the number of rows in L should be
// the same as the lengths of vectors i and w).
vec ccd_kl_update (const mat& L, const vec& u, const vec& w,
		   const vec& x0, double e) {
  unsigned int n = w.n_elem;
  unsigned int m = L.n_cols;
  vec    x  = x0;
  vec    Lx = L * x;
  vec    l(n);
  double h, g, xj, xjnew;
  for (unsigned int j = 0; j < m; j++) {
    l     = L.col(j);
    xj    = x(j);
    compute_grad_ccd_sparse(l,w,Lx,u(j),xj,g,h,e);
    xjnew = project_iterate_ccd(xj,g,h,e);
    Lx   += (xjnew - xj) * l;
    x(j)  = xjnew;
  }
  return x;
}
