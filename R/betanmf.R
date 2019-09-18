# This function decomposes the input matrix X = A*B by nonnegative
# matrix factorization (NMF) based on the beta-divergence criterion
# (negative Poisson log-likelihood) and multiplicative update
# rules. All entries of initial estimates A and B should be
# positive. This is adapted from the MATLAB code by D. Kitamura
# (http://d-kitamura.net).
#
# This function is mainly intended for internal use---specifically, to
# verify implementation of the alternating SQP algorithm.
#
betanmf <- function (X, A, B, numiter = 1000, e = 1e-15) {
  if (!is.matrix(X))
    stop(paste("Input argument \"X\" must be a matrix; see help(matrix)",
               "for more information"))
  for (i in 1:numiter) {

    # Update the loadings ("activations").
    A <- scale.cols(A * ((X / (A %*% B)) %*% t(B)),1/rowSums(B))
    A <- pmax(A,e)

    # Update the factors ("basis vectors").
    B <- B * (t(A) %*% (X / (A %*% B))) / colSums(A)
    B <- pmax(B,e)
  }
  return(list(A = A,B = B))
}
