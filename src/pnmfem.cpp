#include "poismixem.h"

using namespace arma;

// FUNCTION DEFINITIONS
// --------------------
// TO DO: Explain here what this function does, and how to use it.
void pnmfem_update_factor (const mat& X, mat& F, const arma::mat& L,
			   uint j, double numiter) {
  F.col(j) = poismixem(L,X.col(j),F.col(j),numiter);
}

// TO DO: Explain here what this function does, and how to use it.
//
// [[Rcpp::export]]
arma::mat pnmfem_update_factors_rcpp (const arma::mat& X, const arma::mat& F,
				      const arma::mat& L, double numiter) {
  uint m = X.n_cols;
  mat  Fnew = F;
  for (uint j = 0; j < m; j++)
    pnmfem_update_factor(X,Fnew,L,j,numiter);
  return Fnew;
}
