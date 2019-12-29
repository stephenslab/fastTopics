# This function implements the multiplicative update rule for the
# loadings matrix (the "activations"), A, in which the matrix X is
# approximated by the matrix product A*B. Inputs X, A and B should
# not be sparse matrices ("is.matrix" should return TRUE). Input
# argument "e" is a non-negative scalar specifying the minimum value
# of the updated loadings. A positive value of "e" promotes better
# convergence of the multiplicative updates.
betanmf_update_loadings <- function (X, A, B, e = 1e-15) {
  A <- scale.cols(A * tcrossprod(X / (A %*% B),B),1/rowSums(B))
  return(pmax(A,e))
}

# This function implements the multiplicative update rule for the
# factors matrix (the "basis vectors"), B, in which the matrix X is
# approximated by the matrix product A*B. Inputs X, A and B should not
# be sparse matrices ("is.matrix" should return TRUE). Input argument
# "e" is a non-negative scalar specifying the minimum value of the
# updated factors. A positive value of "e" promotes better
# convergence of the multiplicative updates.
betanmf_update_factors <- function (X, A, B, e = 1e-15) {
  B <- B * crossprod(A,X / (A %*% B)) / colSums(A)
  return(pmax(B,e))
}
