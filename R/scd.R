# This function implements a sequential co-ordinate descent (SCD)
# update for the factors matrix (the "basis vectors"), H, in which the
# matrix A is approximated by the matrix product W*H. Inputs A, W and
# H should not be sparse matrices ("is.matrix" should return TRUE).
# Input "numiter" specifies the number of inner-loop iterations to
# perform. Input argument "e" a non-negative scalar included in the
# computations to prevent NaNs due to division by zero.
scd_update_factors <- function (A, W, H, numiter = 1, e = 1e-15)
  scd_update_factors_rcpp(A,t(W),H,numiter,e)

# This function implements a sequential co-ordinate descent (SCD)
# update for the loadings matrix (the "activations"), W, in which the
# matrix A is approximated by the matrix product W*H. Inputs A, W and
# H should not be sparse matrices ("is.matrix" should return TRUE).
# Input "numiter" specifies the number of inner-loop iterations to
# perform. Input argument "e" a non-negative scalar included in the
# computations to prevent NaNs due to division by zero.
scd_update_loadings <- function (A, W, H, numiter = 1, e = 1e-15) {
  A <- t(A)
  W <- t(W)
  return(t(scd_update_factors_rcpp(A,H,W,numiter,e)))
}
