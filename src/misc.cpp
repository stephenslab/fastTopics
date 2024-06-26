#include "misc.h"

using namespace arma;

// FUNCTION DECLARATIONS
// ---------------------
void le_diff (const vec& x, vec& y);

// FUNCTION DEFINITIONS
// --------------------
// Compute, for each row of X, the "least extreme" differences. This
// should output the same result as t(apply(X,1,le.diff)), but faster.
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat le_diff_rcpp (const arma::mat& X) {
  unsigned int n = X.n_rows;
  unsigned int m = X.n_cols;
  mat Y(n,m);
  vec x(m);
  vec y(m);
  for (unsigned int i = 0; i < n; i++) {
    x = trans(X.row(i));
    le_diff(x,y);
    Y.row(i) = trans(y);
  }
  return Y;
}

// This is used to implement x_over_tcrossprod.
// 
// [[Rcpp::export]]
arma::vec x_over_crossprod_rcpp (const arma::vec& i, const arma::vec& j,
				 const arma::vec& x, const arma::mat& A, 
				 const arma::mat& B, double e) {
  unsigned int n = x.n_elem;
  vec y = x;
  for (unsigned int t = 0; t < n; t++)
    y(t) /= (dot(A.col(i(t)),B.col(j(t))) + e);
  return y;
}

// For vector x, return a vector of the same length y containing the
// "least extreme" differences y(i) = x(i) - x(j), in which j is the
// index not equal to i such that abs(x(i) - x(j)) is the smallest
// possible.
void le_diff (const vec& x, vec& y) {
  unsigned int n = x.n_elem;
  if (n == 2) {
    y(0) = x(0) - x(1);
    y(1) = -y(0);
  } else {
    uvec indices = sort_index(x);
    unsigned int i, j, k;
    double a, b;
    i = indices(0);
    j = indices(1);
    y(i) = x(i) - x(j);
    i = indices(n-1);
    j = indices(n-2);
    y(i) = x(i) - x(j);
    for (unsigned int t = 1; t < n-1; t++) {
      i = indices(t-1);
      j = indices(t);
      k = indices(t+1);
      a = x(j) - x(i);
      b = x(k) - x(j);
      if (a <= b)
        y(j) = x(j) - x(i);
      else
        y(j) = x(j) - x(k);
    }
  }
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
//   vec a = nonzeros(A.col(j));
//   unsigned int n = a.n_elem;
//   uvec i(n);
//   getcolnonzeros(A,i,j);
//
void getcolnonzeros (const sp_mat& A, uvec& i, unsigned int j) {
  sp_mat::const_col_iterator ai = A.begin_col(j);
  sp_mat::const_col_iterator an = A.end_col(j);
  for (unsigned int t = 0; ai != an; ++ai, ++t)
    i(t) = ai.row();
}

// Scale each column A[,i] by b[i].
void scalecols (mat& A, const vec& b) {
  rowvec c = trans(b);
  A.each_row() %= c;
}

// Normalize each row of A so that the entries in each row sum to 1.
void normalizerows (mat& A) {
  vec b = conv_to<vec>::from(sum(A,1));
  A.each_col() /= b;
}

// Normalize each column of A so that the entries in each column sum to 1.
void normalizecols (mat& A) {
  rowvec b = sum(A,0);
  A.each_row() /= b;
}

// Scale each row of A so that the largest entry in each row is 1.
void normalizerowsbymax (mat& A) {
  vec b = conv_to<vec>::from(max(A,1));
  A.each_col() /= b;
}
