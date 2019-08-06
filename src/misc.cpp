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

// Copy the ith column of matrix X into previously initialized vector y.
void copycol (const arma::mat& X, int i, arma::vec& y) {
  memcpy(y.memptr(),X.colptr(i),y.n_elem*sizeof(double));
}
