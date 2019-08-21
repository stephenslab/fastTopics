#include "misc.h"

using namespace arma;

// FUNCTION DEFINITIONS
// --------------------
// Return a or b, which ever is smaller.
double min (double a, double b) {
  double y;
  if (a < b)
    y = a;
  else
    y = b;
  return y;
}

// Copy all elements X[i,j] into a previously initialized vector y
// (which should have the same length as i).
void copycolelems (const mat& X, const uvec& i, uint j, vec& y) {
  uint n = i.n_elem;
  for (uint t = 0; t < n; t++)
    y(t) = X(i(t),j);
}

// Copy all elements X[i,j] into a previously initialized vector y
// (which should have the same length as j).
void copyrowelems (const mat& X, uint i, const uvec& j, vec& y) {
  uint m = j.n_elem;
  for (uint t = 0; t < m; t++)
    y(t) = X(i,j(t));
}

// Add b[i] to each row A[i,].
void addtorows (mat& A, const vec& b) {
  uint n = A.n_rows;
  for (uint i = 0; i < n; i++)
    A.row(i) += b(i);
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
