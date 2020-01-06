#include "nnlm.h"

void update (mat& H, const mat & Wt, const mat& A, uint max_iter) {
  // A = W H, solve H
  // No missing in A, Wt = W^T

  uint m = A.n_cols;

  mat WtW;
  vec mu;
  const vec sumW = sum(Wt,1);

  // by columns of H
  for (uint j = 0; j < m; j++)
    scd_kl_update(H.col(j), Wt, A.col(j), sumW, max_iter);
}
