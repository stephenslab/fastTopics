#include "poismix.h"
#include "misc.h"
#include "mixem.h"

using namespace arma;

// INLINE FUNCTION DEFINITIONS
// ---------------------------
// Compute the first-order (g) and second-order (h) partial
// derivatives of the loss function (i.e., the negative Poisson
// log-likelihood) with respect to the jth mixture weight, x[j],
// where: input argument l = L[,j] is the jth column of the n x m data
// matrix, L; input w is a vector of length n containing a count or
// "pseudocount" associated with each row of L; x is the initial
// estimate of the mixture weights; and Lx = L*x. Input argument u is
// a vector storing an intermediate result of the same size as w.
//
// This is used by scd_kl_update.
//
inline void compute_grad_kl (const vec& l, const vec& w, const vec& Lx,
			     double x, double& g, double& h, vec& r,
			     double e) {
  r = l/(Lx + e);
  h = dot(square(r),w);
  g = dot(r,w) - sum(l) + h*x;
}

// This is the same as compute_grad_kl, except that the indices of the
// nonzero counts are supplied in vector i, and w must contain (i and
// w should be the same length).
inline void compute_grad_kl_sparse (const vec& l, const vec& w, const uvec& i,
				    double x, const vec& Lx, double& g,
				    double& h, vec& r, double e) {
  r = l/(Lx + e);
  h = dot_square_sparse_b(r,w,i);
  g = dot_sparse_b(r,w,i) - sum(l) + h*x;
}

// Given the first-order (g) and second-order (h) partial derivatives
// of the loss function (negative Poisson log-likelihood), compute the
// updated value of the mixture weight x[j] suitably projected so that
// it is non-negative.
//
// This is used by scd_kl_update.
//
inline double project_iterate_kl (double g, double h, double e) {
  double x = g/(h + e);
  x = maximum(x,0);
  return x;
}

// FUNCTION DEFINITIONS
// --------------------
// This is mainly used to test the first variant of the poismixem C++
// function.
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec poismixem_rcpp (const arma::mat& L, const arma::vec& w,
			  const arma::vec& x0, uint numiter) {
  return poismixem(L,w,x0,numiter);
}

// This is mainly used to test the first variant of the poismixsqp C++
// interface.
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec poismixsqp_rcpp (const arma::mat& L, const arma::vec& w,
			   const arma::vec& x0, uint numiter,
			   const Rcpp::List& control) {
  mixsqp_control_params ctrl = get_mixsqp_control_params(control);
  return poismixsqp(L,w,x0,numiter,ctrl);
}

// This is mainly used to test the second variant of the poismixem C++
// interface.
//
// [[Rcpp::export]]
arma::vec poismixem2_rcpp (const arma::mat& L1, const arma::vec& w,
			   const arma::vec& u, const arma::vec& x0,
			   uint numiter) {
  vec x = x0;
  mat P = L1;
  poismixem(L1,u,w,x,P,numiter);
  return x;
}

// This is mainly used to test the second variant of the poismixsqp C++
// interface.
//
// [[Rcpp::export]]
arma::vec poismixsqp2_rcpp (const arma::mat& L1, const arma::vec& w,
			    const arma::vec& u, const arma::vec& x0,
			    uint numiter, const Rcpp::List& control) {
  mixsqp_control_params ctrl = get_mixsqp_control_params(control);
  uint m = L1.n_cols;
  vec  x = x0;
  mat  Z = L1;
  mat  H(m,m);
  poismixsqp(L1,u,w,x,Z,H,numiter,ctrl);
  return x;
}

// This is mainly used to test the third variant of the poismixem C++
// interface.
//
// [[Rcpp::export]]
arma::vec poismixem3_rcpp (const arma::mat& L1, const arma::vec& w,
			   const arma::vec& u, const arma::uvec& i,
			   const arma::vec& x0, uint numiter) {
  vec x = x0;
  poismixem(L1,u,w,i,x,numiter);
  return x;
}

// This is mainly used to test the third variant of the poismixsqp C++
// interface.
//
// [[Rcpp::export]]
arma::vec poismixsqp3_rcpp (const arma::mat& L1, const arma::vec& w,
			    const arma::vec& u, const arma::uvec& i,
			    const arma::vec& x0, uint numiter,
			    const Rcpp::List& control) {
  mixsqp_control_params ctrl = get_mixsqp_control_params(control);
  uint m = L1.n_cols;
  vec  x = x0;
  mat  H(m,m);
  poismixsqp(L1,u,w,i,x,H,numiter,ctrl);
  return x;
}

// This is mainly used to test C++ function scd_kl_update.
//
// [[Rcpp::export]]
arma::vec scd_kl_update_rcpp (const arma::mat& L, const arma::vec& w, 
			      const arma::vec& x0, uint numiter, double e) {
  return scd_kl_update(L,w,x0,numiter,e);
}

// This is mainly used to test C++ function scd_kl_update_sparse.
//
// [[Rcpp::export]]
arma::vec scd_kl_update_sparse_rcpp (const arma::mat& L, const arma::vec& w, 
				     const arma::uvec& i, const arma::vec& x0,
				     uint numiter, double e) {
  return scd_kl_update_sparse(L,w,i,x0,numiter,e);
}

// NOTES:
//  + L is is a m x n matrix.
//  + w is a vector of length n.
//  + x is a vector of length m; it will be updated.
//  + Lx = L' * x, and will be updated.
//
// This is mainly used to test C++ function ccd_kl_update.
//
// [[Rcpp::export]]
void ccd_kl_update_rcpp (const Rcpp::NumericMatrix& L,
			 const Rcpp::NumericVector& w,
			 Rcpp::NumericVector& Lx,
			 Rcpp::NumericVector& x,
			 double e) {
  ccd_kl_update(L.ncol(),L.nrow(),x.begin(),Lx.begin(),w.begin(),L.begin(),e);
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
vec poismixem (const mat& L, const vec& w, const vec& x0, uint numiter) {
  mat L1 = L;
  mat P  = L;
  vec u  = sum(L,0);
  vec x  = x0;
  normalizecols(L1);
  poismixem(L1,u,w,x,P,numiter);
  return x;
}

// Compute a maximum-likelihood estimate (MLE) of the mixture weights
// in a Poisson mixture model by iterating the mix-SQP updates for a
// fixed number of iterations.
//
// Input argument L is an n x m matrix with non-negative entries;
// input w is a vector of length n containing a count or "pseudocount"
// associated with each row of L; input argument x0 is the initial
// estimate of the mixture weights; and input "numiter" specifies the
// number of mix-SQP updates to perform.
//
// The return value is a vector of length m containing the updated
// mixture weights.
vec poismixsqp (const mat& L, const vec& w, const vec& x0, uint numiter,
		const mixsqp_control_params& control) {
  int m  = L.n_cols;
  mat L1 = L;
  mat Z  = L;
  vec x  = x0;
  vec u  = sum(L,0);
  mat H(m,m);
  normalizecols(L1);
  poismixsqp(L1,u,w,x,Z,H,numiter,control);
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
		uint numiter) {
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

// Use this variant of poismixsqp if you plan on calling poismixsqp
// multiple times with the same matrix, L. In this call, input u should
// contain the column sums, u = colSums(L), in which L is the matrix
// prior to normalization, and matrix L1 is the normalized version of
// L in which each column sums to 1; that is, L1 = normalize.cols(L).
// Input Z should be a matrix of the same size as L1, and H should be
// an m x m matrix, where m is the number of columns in L.
//
// Note that in this variant of poismixem, L1 and w do not need to
// contain all the data; once u is correctly calculated, any rows of
// L1 associated with zero weights have no effect, so only the vector
// of nonzero weights w, and the rows of L1 associated with those
// weights, need to be supplied.
void poismixsqp (const mat& L1, const vec& u, const vec& w, vec& x, mat& Z,
		 mat& H, uint numiter, const mixsqp_control_params& control) {
  uvec i = find(w > 0);
  if (i.n_elem == 1) {

    // Treat the special case when only one of the counts is nonzero.
    poismix_one_nonzero(L1,u,w,i(0),x);
  } else {
    vec    objective(numiter);
    double s = sum(w);

    // Recover the mixture proportions of the multinomial mixture model
    // from the mixture weights of the Poisson mixture model.
    x %= u;
    x /= s;

    // Perform one or more SQP updates for the multinomial mixture model.
    mixsqp(L1,w,x,Z,H,numiter,control,objective);

    // Recover the mixture weights of the Poisson mixture model from the
    // mixture weights of the multinomial mixture model.
    x *= s;
    x /= u;
  }
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
		vec& x, uint numiter) {
  uint n = i.n_elem;
  uint m = x.n_elem;
  mat  P(n,m);
  poismixem(L1.rows(i),u,w,x,P,numiter);
}

// The third poismixsqp interface is similar to the second, except
// that the indices of the nonzero counts are supplied in vector i,
// and w must contain only the nonzero counts (i and w should be the
// same length).
//
// To illustrate, in this example x1 and x2 should be the same after
// running this R code:
//
//   i  <- which(w > 0)
//   x1 <- poismixsqp2_rcpp(L1,w,u,x0,numiter,control)
//   x2 <- poismixsqp3_rcpp(L1,w[i],u,i-1,x0,numiter,control)
//
void poismixsqp (const mat& L1, const vec& u, const vec& w, const uvec& i,
		 vec& x, mat& H, uint numiter,
		 const mixsqp_control_params& control) {
  uint n = i.n_elem;
  uint m = x.n_elem;
  mat  Z(n,m);
  vec  objective(numiter);
  poismixsqp(L1.rows(i),u,w,x,Z,H,numiter,control);
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
		   uint numiter, double e) {
  uint   n  = L.n_rows;
  uint   m  = L.n_cols;
  vec    x  = x0;
  vec    Lx = L * x;
  vec    l(n);
  vec    r(n);
  double h, g, xj, xjnew;
  for (uint iter = 0; iter < numiter; iter++)
    for (uint j = 0; j < m; j++) {
      l     = L.col(j);
      xj    = x(j);
      compute_grad_kl(l,w,Lx,xj,g,h,r,e);
      xjnew = project_iterate_kl(g,h,e);
      Lx   += (xjnew - xj) * l;
      x(j)  = xjnew;
    }
  return x;
}

// This is the same as scd_kl_update, except that the indices of the
// nonzero counts are supplied in vector i, and w must contain (i and
// w should be the same length).
vec scd_kl_update_sparse (const mat& L, const vec& w, const uvec& i,
			  const vec& x0, uint numiter, double e) {
  uint   n  = L.n_rows;
  uint   m  = L.n_cols;
  vec    x  = x0;
  vec    Lx = L * x;
  vec    l(n);
  vec    r(n);
  double h, g, xj, xjnew;
  for (uint iter = 0; iter < numiter; iter++)
    for (uint j = 0; j < m; j++) {
      l     = L.col(j);
      xj    = x(j);
      compute_grad_kl_sparse(l,w,i,xj,Lx,g,h,r,e);
      xjnew = project_iterate_kl(g,h,e);
      Lx   += (xjnew - xj) * l;
      x(j)  = xjnew;
    }
  return x;
}

// Implements the core part of the cyclic co-ordinate descent (CCD)
// updates.
void ccd_kl_update (uint m, uint k, double* wt, double* wht,
		    const double* vt, const double* H, double e) {
  double d, g, h, t, w0, w1;
  for (uint i = 0; i < k; i++) {
    g = 0;
    h = 0;
    for (uint j = 0, hi = i; j < m; j++, hi += k) {
      t  = vt[j]/(wht[j] + e);
      g += H[hi]*(1 - t);
      h += H[hi]*H[hi]*t/(wht[j] + e);
    }
    w0    = wt[i];
    w1    = wt[i] - g/h + e;
    w1    = maximum(w1,e);
    d     = w1 - w0;
    wt[i] = w1;
    for (uint j = 0; j < m; j++)
      wht[j] += d * H[j*k + i];
  }
}

// Find the maximum-likelihood estimate (MLE) for the special case
// when only one of the counts is positive. Input argument w should be
// a vector of counts (in which exactly one of these counts is
// positive), and i should be the index of the nonzero count. See
// above for an explanation of inputs L1 and u.
void poismix_one_nonzero (const mat& L1, const vec& u, const vec& w, 
			  uint i, vec& x) {
  mixture_one_nonzero(L1,i,x);
  uint j = index_max(x);
  x(j)   = w(i)/u(j);
}
