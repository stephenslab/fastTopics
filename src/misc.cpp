#include "misc.h"

using namespace arma;

// FUNCTION DEFINITIONS
// --------------------

// This is used to implement x_over_tcrossprod.
// 
// [[Rcpp::export]]
arma::vec x_over_crossprod_rcpp (const arma::vec& i, const arma::vec& j,
				 const arma::vec& x, const arma::mat& A, 
				 const arma::mat& B, double e) {
  uint n = x.n_elem;
  vec  y = x;
  for (uint t = 0; t < n; t++)
    y(t) /= (dot(A.col(i(t)),B.col(j(t))) + e);
  return y;
}

// Return the row indices of the nonzeros in the jth column of sparse
// matrix A. This is the same as
//
//  i = find(A.col(j))
//
// but this code does not compile in some versions of gcc, so I
// re-implemented this code here. Vector i must already been
// initialized with the proper length, e.g., by doing
//
//   vec  a = nonzeros(A.col(j));
//   uint n = a.n_elem;
//   uvec i(n);
//   getcolnonzeros(A,i,j);
//
void getcolnonzeros (const sp_mat& A, uvec& i, uint j) {
  sp_mat::const_col_iterator ai = A.begin_col(j);
  sp_mat::const_col_iterator an = A.end_col(j);
  for (uint t = 0; ai != an; ++ai, ++t)
    i(t) = ai.row();
}

// Scale each column A[,i] by b[i].
void scalecols (mat& A, const vec& b) {
  vec c = b;
  A.each_row() %= c;
}

// Normalize each row of A so that the entries in each row sum to 1.
void normalizerows (mat& A) {
  vec b = sum(A,1);
  A.each_col() /= b;
}

// Normalize each row of A so that the entries in each column sum to 1.
void normalizecols (mat& A) {
  vec b = sum(A,0);
  A.each_row() /= b;
}

// Scale each row of A so that the largest entry in each row is 1.
void normalizerowsbymax (mat& A) {
  vec b = max(A,1);
  A.each_col() /= b;
}

