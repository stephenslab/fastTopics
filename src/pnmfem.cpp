#include "poismixem.h"

using namespace arma;

// FUNCTION DEFINITIONS
// --------------------
// Update the jth column of the factors matrix, F. See function
// pnmfem_update_factors_rcpp for further details.
//
// Input P should be a matrix of the same dimension as L. It is used
// to store a temporary result. Also note that L is modified in the
// process of performing the update, and should not be reused.
inline void pnmfem_update_factor (const mat& X, mat& F, mat& L, mat& P, 
				  uint j, double numiter) {
  // F.col(j) = poismixem(L,X.col(j),F.col(j),P,numiter);
}

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
  uint m = X.n_cols;
  mat  Fnew = F;
  mat  L1   = L;
  mat  P    = L;
  for (uint j = 0; j < m; j++) {
    L1 = L;
    pnmfem_update_factor(X,Fnew,L1,P,j,numiter);
  }
  return Fnew;
}
