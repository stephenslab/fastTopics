#include <RcppParallel.h>
#include "misc.h"
#include "poismix.h"

using namespace arma;

// INLINE FUNCTION DEFINITIONS
// ---------------------------
// Perform one or several EM updates for a single column of the k x m
// factors matrix, F, when X is a dense matrix.
inline vec pnmfem_update_factor (const mat& X, const mat& F, const mat& L1,
				 const vec& u, mat& P, unsigned int j, 
				 unsigned int numiter) {
  colvec temp = F.col(j);
  vec f = F.col(j);
  poismixem(L1,u,X.col(j),f,P,numiter);
  return f;
}

// Perform one or several EM updates for a single column of the k x m
// factors matrix, F, when X is a sparse matrix.
inline vec pnmfem_update_factor_sparse (const sp_mat& X, const mat& F,
					const mat& L1, const vec& u, 
					unsigned int j, unsigned int numiter) {
  vec          x = nonzeros(X.col(j));
  unsigned int n = x.n_elem;
  vec          f = F.col(j);
  uvec         i(n);
  getcolnonzeros(X,i,j);
  poismixem(L1,u,x,i,f,numiter);
  return f;
}

// CLASS DEFINITIONS
// -----------------
// This class is used to implement multithreaded computation of the
// factor updates in pnmfem_update_factors_parallel_rcpp.
struct pnmfem_factor_updater : public RcppParallel::Worker {
  const mat&   X;
  const mat&   F;
  mat          L1;
  vec          u;
  mat&         Fnew;
  const vec&   j;
  unsigned int numiter;

  // This is used to create a pnmfem_factor_updater object.
  pnmfem_factor_updater (const mat& X, const mat& F, const mat& L, 
			 mat& Fnew, const vec& j, unsigned int numiter) :
    X(X), F(F), L1(L), u(L.n_cols), Fnew(Fnew), j(j), numiter(numiter) { 
    u = trans(sum(L,0));
    normalizecols(L1);
  };

  // This function updates the factors for a given range of columns.
  void operator() (std::size_t begin, std::size_t end) {
    mat P = L1;
    for (unsigned int i = begin; i < end; i++)
      Fnew.col(j(i)) = pnmfem_update_factor(X,F,L1,u,P,j(i),numiter);
  }
};

// This class is used to implement multithreaded computation of the
// factor updates in pnmfem_update_factors_sparse_parallel_rcpp.
struct pnmfem_factor_updater_sparse : public RcppParallel::Worker {
  const sp_mat& X;
  const mat&    F;
  mat           L1;
  vec           u;
  mat&          Fnew;
  const vec&    j;
  unsigned int  numiter;

  // This is used to create a pnmfem_factor_updater object.
  pnmfem_factor_updater_sparse (const sp_mat& X, const mat& F, const mat& L, 
				mat& Fnew, const vec& j, 
				unsigned int numiter) :
    X(X), F(F), L1(L), u(L.n_cols), Fnew(Fnew), j(j), numiter(numiter) {
    u = trans(sum(L,0));
    normalizecols(L1);
  };

  // This function updates the factors for a given range of columns.
  void operator() (std::size_t begin, std::size_t end) {
    for (unsigned int i = begin; i < end; i++)
      Fnew.col(j(i)) = pnmfem_update_factor_sparse(X,F,L1,u,j(i),numiter);
  }
};

// FUNCTION DEFINITIONS
// --------------------
// Perform one or several EM updates for the factors matrix, F, in
// which the matrix X is approximated by L*F. Input "k" specifies
// which columns of F are updated. Input "numiter" specifies the
// number of EM updates to perform.
//
// Note that, unlike most other functions implemented in this package,
// the factors matrix F should be a k x m matrix, where k is the
// number factors, or "topics". This is done for ease of
// implementation, and for speed, because the computation is performed
// on the m columns of F, and the matrix is stored columnwise.
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat pnmfem_update_factors_rcpp (const arma::mat& X, const arma::mat& F,
				      const arma::mat& L, const arma::vec& j,
				      double numiter) {
  unsigned int n = j.n_elem;
  vec  u    = trans(sum(L,0));
  mat  L1   = L;
  mat  P    = L;
  mat  Fnew = F;
  normalizecols(L1);
  for (unsigned int i = 0; i < n; i++)
    Fnew.col(j(i)) = pnmfem_update_factor(X,F,L1,u,P,j(i),numiter);
  return Fnew;
}

// This does the same thing as pnmfem_update_factors_rcpp, except that
// X is a sparse matrix. See pnmfem_update_factors_rcpp for details.
//
// [[Rcpp::export]]
arma::mat pnmfem_update_factors_sparse_rcpp (const arma::sp_mat& X,
					     const arma::mat& F,
					     const arma::mat& L,
					     const arma::vec& j,
					     double numiter) {
  unsigned int n = j.n_elem;
  vec  u    = trans(sum(L,0));
  mat  L1   = L;
  mat  Fnew = F;
  normalizecols(L1);
  for (unsigned int i = 0; i < n; i++)
    Fnew.col(j(i)) = pnmfem_update_factor_sparse(X,F,L1,u,j(i),numiter);
  return Fnew;
}

// This does the same thing as pnmfem_update_factors_rcpp, except that
// Intel Threading Building Blocks (TBB) are used to update the
// factors in parallel.
//
// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::export]]
arma::mat pnmfem_update_factors_parallel_rcpp (const arma::mat& X,
					       const arma::mat& F,
					       const arma::mat& L,
					       const arma::vec& j,
					       double numiter) {
  mat Fnew = F;
  pnmfem_factor_updater worker(X,F,L,Fnew,j,numiter);
  parallelFor(0,j.n_elem,worker);
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
						      const arma::vec& j,
						      double numiter) {
  mat Fnew = F;
  pnmfem_factor_updater_sparse worker(X,F,L,Fnew,j,numiter);
  parallelFor(0,j.n_elem,worker);
  return Fnew;
}
