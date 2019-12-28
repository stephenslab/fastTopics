#include "misc.h"
#include "poismixem.h"

using namespace arma;

// FUNCTION DEFINITIONS
// --------------------
// Perform an EM update for the factors matrix, F, in which the matrix
// X is approximated by L*F. Input "numiter" specifies the number of
// EM updates to perform.
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
  mat  Fnew = F;
  mat  L1   = L;
  mat  P    = L;
  vec  u    = sum(L,0);
  vec  f(m);
  normalizecols(L1);
  for (uint j = 0; j < m; j++) {
    f = Fnew.col(j);
    poismixem(L1,u,X.col(j),f,P,numiter);
    Fnew.col(j) = f;
  }
  return Fnew;
}

// TO DO: Explain here what this function does, and how to use it.
//
// [[Rcpp::export]]
arma::mat pnmfem_update_factors_sparse_rcpp (const arma::sp_mat& X,
					     const arma::mat& F,
					     const arma::mat& L,
					     double numiter) {
  mat Fnew = F;
  // TO DO.
  return Fnew;
}
