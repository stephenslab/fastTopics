#include <RcppParallel.h>
#include "misc.h"
#include "poismixem.h"

using namespace arma;

// INLINE FUNCTION DEFINITIONS
// ---------------------------
// Perform one or several EM updates for a single column of the k x m
// factors matrix, F.
inline vec pnmfem_update_factor (const mat& X, const mat& F, const mat& L1,
				 const vec& u, mat& P, uint j, uint numiter) {
  vec f = F.col(j);
  poismixem(L1,u,X.col(j),f,P,numiter);
  return f;
}

// Perform one or several EM updates for a single column of the k x m
// factors matrix, F.
inline vec pnmfem_update_factor_sparse (const sp_mat& X, const mat& F,
					const mat& L1, const vec& u, uint j,
					uint numiter) {
  uint k = F.n_rows;
  vec  f = F.col(j);
  vec  x = nonzeros(X.col(j));
  uint n = x.n_elem;
  uvec i(n);
  mat  P(n,k);
  getcolnonzeros(X,i,j);
  poismixem(L1.rows(i),u,x,f,P,numiter);
  return f;
}

// CLASS DEFINITIONS
// -----------------
// This class is used to implement multithreaded computation of the
// factor updates in pnmfem_update_factors_parallel_rcpp.
//
// [[Rcpp::depends(RcppParallel)]]
struct pnmfem_factor_updater : public RcppParallel::Worker {
  const mat& X;
  const mat& F;
  const mat& L1;
  const vec& u;
  mat&       Fnew;
  uint       numiter;

  // This is used to create a pnmfem_factor_updater object.
  pnmfem_factor_updater (const arma::mat& X, const arma::mat& F,
			 const arma::mat& L1, const vec& u, mat& Fnew,
			 uint numiter) :
    X(X), F(F), L1(L1), u(u), Fnew(Fnew), numiter(numiter) { };

  // This function updates the factors for a given range of columns.
  void operator() (std::size_t begin, std::size_t end) {
    mat P = L1;
    for (uint j = begin; j < end; j++)
      Fnew.col(j) = pnmfem_update_factor(X,F,L1,u,P,j,numiter);
  }
};

// This class is used to implement multithreaded computation of the
// factor updates in pnmfem_update_factors_sparse_parallel_rcpp.
//
// [[Rcpp::depends(RcppParallel)]]
struct pnmfem_factor_updater_sparse : public RcppParallel::Worker {
  const sp_mat& X;
  const mat&    F;
  const mat&    L1;
  const vec&    u;
  mat&          Fnew;
  uint          numiter;

  // This is used to create a pnmfem_factor_updater object.
  pnmfem_factor_updater_sparse (const arma::sp_mat& X, const arma::mat& F,
				const arma::mat& L1, const vec& u, mat& Fnew,
				uint numiter) :
    X(X), F(F), L1(L1), u(u), Fnew(Fnew), numiter(numiter) { };

  // This function updates the factors for a given range of columns.
  void operator() (std::size_t begin, std::size_t end) {
    for (uint j = begin; j < end; j++)
      Fnew.col(j) = pnmfem_update_factor_sparse(X,F,L1,u,j,numiter);
  }
};

// FUNCTION DEFINITIONS
// --------------------
// Perform one or several EM updates for the factors matrix, F, in
// which the matrix X is approximated by L*F. Input "numiter"
// specifies the number of EM updates to perform.
//
// Note that, unlike most other functions implemented in this package,
// the factors matrix F should be a k x m matrix, where k is the
// number factors, or "topics". This is done for ease of
// implementation, and for speed, because the computation is performed
// on the m columns of F, and the matrix is stored columnwise.
//
// [[Rcpp::export]]
arma::mat pnmfem_update_factors_rcpp (const arma::mat& X, const arma::mat& F,
				      const arma::mat& L, double numiter) {
  uint m    = X.n_cols;
  vec  u    = sum(L,0);
  mat  L1   = L;
  mat  P    = L;
  mat  Fnew = F;
  normalizecols(L1);
  for (uint j = 0; j < m; j++)
    Fnew.col(j) = pnmfem_update_factor(X,F,L1,u,P,j,numiter);
  return Fnew;
}

// This does the same thing as pnmfem_update_factors_rcpp, except that
// X is a sparse matrix. See pnmfem_update_factors_rcpp for details.
//
// [[Rcpp::export]]
arma::mat pnmfem_update_factors_sparse_rcpp (const arma::sp_mat& X,
					     const arma::mat& F,
					     const arma::mat& L,
					     double numiter) {
  uint m    = X.n_cols;
  vec  u    = sum(L,0);
  mat  L1   = L;
  mat  Fnew = F;
  normalizecols(L1);
  for (uint j = 0; j < m; j++)
    Fnew.col(j) = pnmfem_update_factor_sparse(X,F,L1,u,j,numiter);
  return Fnew;
}

// This does the same thing as pnmfem_update_factors_rcpp, except that
// Intel Threading Building Blocks (TBB) are used to update the
// factors in parallel.
//
// [[Rcpp::export]]
arma::mat pnmfem_update_factors_parallel_rcpp (const arma::mat& X,
					       const arma::mat& F,
					       const arma::mat& L,
					       double numiter) {
  uint m    = X.n_cols;
  vec  u    = sum(L,0);
  mat  L1   = L;
  mat  Fnew = F;
  normalizecols(L1);
  pnmfem_factor_updater worker(X,F,L1,u,Fnew,numiter);
  parallelFor(0,m,worker);
  return Fnew;
}

// This does the same thing as pnmfem_update_factors_sparse_rcpp,
// except that Intel Threading Building Blocks (TBB) are used to
// update the factors in parallel.
//
// [[Rcpp::export]]
arma::mat pnmfem_update_factors_sparse_parallel_rcpp (const arma::sp_mat& X,
						      const arma::mat& F,
						      const arma::mat& L,
						      double numiter) {
  uint m    = X.n_cols;
  vec  u    = sum(L,0);
  mat  L1   = L;
  mat  Fnew = F;
  normalizecols(L1);
  pnmfem_factor_updater_sparse worker(X,F,L1,u,Fnew,numiter);
  parallelFor(0,m,worker);
  return Fnew;
}
