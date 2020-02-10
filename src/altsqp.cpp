#include <RcppParallel.h>
#include "misc.h"
#include "poismix.h"

using namespace arma;

// INLINE FUNCTION DEFINITIONS
// ---------------------------
// Perform one or more alt-SQP updates for a single column of the k x m
// factors matrix, F, when X is a dense matrix.
inline vec altsqp_update_factor (const mat& X, const mat& F, const mat& L1,
				 const vec& u, mat& Z, mat& H, unsigned int j,
				 unsigned int numiter,
				 const mixsqp_control_params& control) {
  vec f = F.col(j);
  poismixsqp(L1,u,X.col(j),f,Z,H,numiter,control);
  return f;
}

// Perform one or more alt-SQP updates for a single column of the k x m
// factors matrix, F, when X is a sparse matrix.
inline vec altsqp_update_factor_sparse (const sp_mat& X, const mat& F,
					const mat& L1, const vec& u, mat& H,
					unsigned int j, unsigned int numiter,
					const mixsqp_control_params& control) {
  vec          x = nonzeros(X.col(j));
  unsigned int n = x.n_elem;
  vec          f = F.col(j);
  uvec         i(n);
  getcolnonzeros(X,i,j);
  poismixsqp(L1,u,x,i,f,H,numiter,control);
  return f;
}

// CLASS DEFINITIONS
// -----------------
// This class is used to implement multithreaded computation of the
// factor updates in altsqp_update_factors_parallel_rcpp.
struct altsqp_factor_updater : public RcppParallel::Worker {
  const mat&   X;
  const mat&   F;
  mat          L1;
  vec          u;
  mat&         Fnew;
  unsigned int numiter;
  const mixsqp_control_params& control;

  // This is used to create a altsqp_factor_updater object.
  altsqp_factor_updater (const mat& X, const mat& F, const mat& L, 
			 mat& Fnew, unsigned int numiter, 
			 const mixsqp_control_params& control) :
    X(X), F(F), L1(L), u(L.n_cols), Fnew(Fnew), numiter(numiter), 
    control(control) { 
    u = sum(L,0);
    normalizecols(L1);
  };

  // This function updates the factors for a given range of columns.
  void operator() (std::size_t begin, std::size_t end) {
    unsigned int m = X.n_cols;
    mat Z = L1;
    mat H(m,m);
    for (unsigned int j = begin; j < end; j++)
      Fnew.col(j) = altsqp_update_factor(X,F,L1,u,Z,H,j,numiter,control);
  }
};

// This class is used to implement multithreaded computation of the
// factor updates in altsqp_update_factors_sparse_parallel_rcpp.
struct altsqp_factor_updater_sparse : public RcppParallel::Worker {
  const sp_mat& X;
  const mat&    F;
  mat           L1;
  vec           u;
  mat&          Fnew;
  unsigned int  numiter;
  const mixsqp_control_params& control;

  // This is used to create a altsqp_factor_updater object.
  altsqp_factor_updater_sparse (const sp_mat& X, const mat& F, const mat& L, 
				mat& Fnew, unsigned int numiter,
				const mixsqp_control_params& control) :
    X(X), F(F), L1(L), u(L.n_cols), Fnew(Fnew), numiter(numiter), 
    control(control) { 
    u = sum(L,0);
    normalizecols(L1);
  };

  // This function updates the factors for a given range of columns.
  void operator() (std::size_t begin, std::size_t end) {
    unsigned int m = X.n_cols;
    mat H(m,m);
    for (unsigned int j = begin; j < end; j++)
      Fnew.col(j) = altsqp_update_factor_sparse(X,F,L1,u,H,j,numiter,control);
  }
};

// FUNCTION DEFINITIONS
// --------------------
// Perform one or more alternating SQP ("alt-SQP") updates for the
// factors matrix, F, in which the matrix X is approximated by
// L*F. Input "numiter" specifies the number of EM updates to perform.
//
// Note that, unlike most other functions implemented in this package,
// the factors matrix F should be a k x m matrix, where k is the
// number factors, or "topics". This is done for ease of
// implementation, and for speed, because the computation is performed
// on the m columns of F, and the matrix is stored columnwise.
// 
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat altsqp_update_factors_rcpp (const arma::mat& X, const arma::mat& F,
				      const arma::mat& L, double numiter,
				      const Rcpp::List& control) {
  mixsqp_control_params ctrl = get_mixsqp_control_params(control);
  unsigned int m = X.n_cols;
  vec  u    = sum(L,0);
  mat  L1   = L;
  mat  Z    = L;
  mat  Fnew = F;
  mat  H(m,m);
  normalizecols(L1);
  for (unsigned int j = 0; j < m; j++)
    Fnew.col(j) = altsqp_update_factor(X,F,L1,u,Z,H,j,numiter,ctrl);
  return Fnew;
}

// This does the same thing as altsqp_update_factors_rcpp, except that
// X is a sparse matrix. See altsqp_update_factors_rcpp for details.
//
// [[Rcpp::export]]
arma::mat altsqp_update_factors_sparse_rcpp (const arma::sp_mat& X,
					     const arma::mat& F,
					     const arma::mat& L,
					     double numiter,
					     const Rcpp::List& control) {
  mixsqp_control_params ctrl = get_mixsqp_control_params(control);
  unsigned int m = X.n_cols;
  vec  u    = sum(L,0);
  mat  L1   = L;
  mat  Z    = L;
  mat  Fnew = F;
  mat  H(m,m);
  normalizecols(L1);
  for (unsigned int j = 0; j < m; j++)
    Fnew.col(j) = altsqp_update_factor_sparse(X,F,L1,u,H,j,numiter,ctrl);
  return Fnew;
}

// This does the same thing as altsqp_update_factors_rcpp, except that
// Intel Threading Building Blocks (TBB) are used to update the
// factors in parallel.
//
// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::export]]
arma::mat altsqp_update_factors_parallel_rcpp (const arma::mat& X,
					       const arma::mat& F,
					       const arma::mat& L,
					       double numiter,
					       const Rcpp::List& control) {
  mixsqp_control_params ctrl = get_mixsqp_control_params(control);
  mat Fnew = F;
  altsqp_factor_updater worker(X,F,L,Fnew,numiter,ctrl);
  parallelFor(0,X.n_cols,worker);
  return Fnew;
}

// This does the same thing as altsqp_update_factors_sparse_rcpp,
// except that Intel Threading Building Blocks (TBB) are used to
// update the factors in parallel.
//
// [[Rcpp::export]]
arma::mat altsqp_update_factors_sparse_parallel_rcpp (const arma::sp_mat& X,
	    const arma::mat& F, const arma::mat& L, double numiter,
	    const Rcpp::List& control) {
  mixsqp_control_params ctrl = get_mixsqp_control_params(control);
  mat Fnew = F;
  altsqp_factor_updater_sparse worker(X,F,L,Fnew,numiter,ctrl);
  parallelFor(0,X.n_cols,worker);
  return Fnew;
}
