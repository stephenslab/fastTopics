# This function implements the cyclic co-ordinate descent (CCD) update
# for the factors matrix (the "basic vectors"), B, in which the matrix
# X is approximated by the matrix product A*B. Inputs X, A and B
# should not be sparse matrices ("is.matrix" should return
# TRUE). Input argument "e" is a non-negative scalar specifying the
# minimum value of the updated factors.
#
# Note that the call to ccd_update_factors_rcpp modifies inputs B1 and
# AB.
ccd_update_factors <- function (X, A, B, e = 1e-15) {
  A  <- t(A)
  AB <- crossprod(A,B)
  B1 <- B + 0
  ccd_update_factors_rcpp(X,A,B1,AB,e)
  return(B1)
}

# This function implements the cyclic co-ordinate descent (CCD) update
# for the loadings matrix (the "activations"), A, in which the matrix
# X is approximated by the matrix product A*B. Inputs X, A and B
# should not be sparse matrices ("is.matrix" should return
# TRUE). Input argument "e" is a non-negative scalar specifying the
# minimum value of the updated factors.
#
# Note that the call to ccd_update_factors_rcpp modifies inputs A and
# AB.
ccd_update_loadings <- function (X, A, B, e = 1e-15) {
  X  <- t(X)
  A  <- t(A)
  AB <- crossprod(B,A)
  ccd_update_factors_rcpp(X,B,A,AB,e)
  return(t(A))
}
