# This function implements the multiplicative update rule for the
# loadings matrix (the "activations"), A, in which the matrix X is
# approximated by the matrix product A*B. Inputs X, A and B should
# not be sparse matrices ("is.matrix" should return TRUE).
betanmf_update_loadings <- function (X, A, B)
  scale.cols(A * tcrossprod(X / (A %*% B),B),1/rowSums(B))

# This function implements the multiplicative update rule for the
# factors matrix (the "basis vectors"), B, in which the matrix X is
# approximated by the matrix product A*B. Inputs X, A and B should not
# be sparse matrices ("is.matrix" should return TRUE).
betanmf_update_factors <- function (X, A, B)
  B * crossprod(A,X / (A %*% B)) / colSums(A)
