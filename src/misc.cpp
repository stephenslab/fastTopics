#include "misc.h"

using namespace arma;

// FUNCTION DEFINITIONS
// --------------------
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
  uint n = A.n_cols;
  for (uint i = 0; i < n; i++)
    A.col(i) *= b(i);
}

// Normalize each row of A so that the entries in each row sum to 1.
void normalizerows (mat& A) {
  uint n = A.n_rows;
  vec  b = sum(A,1);
  for (uint i = 0; i < n; i++) 
    A.row(i) /= b(i);
}
