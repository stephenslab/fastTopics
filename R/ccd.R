# This function implements the cyclic co-ordinate descent (CCD) update
# for the factors matrix (the "basis vectors"), B, in which the matrix
# V is approximated by the matrix product W*H. Inputs W and H should
# be dense matrices ("is.matrix" should return TRUE). Input argument
# "e" is a non-negative scalar specifying the minimum value of the
# updated factors.
#
#' @importFrom Rcpp evalCpp
#'
ccd_update_factors <- function (V, W, H, e = 1e-15) {
  H <- ccd_update_factors_rcpp(V,W,H,e)
  return(H)
}

# This function implements the cyclic co-ordinate descent (CCD) update
# for the loadings matrix (the "activations"), V, in which the matrix
# V is approximated by the matrix product W*H. Inputs W and H should
# not be sparse matrices ("is.matrix" should return TRUE). Input
# argument "e" is a non-negative scalar specifying the minimum value
# of the updated factors.
#
#' @importFrom Rcpp evalCpp
#'
ccd_update_loadings <- function (V, W, H, e = 1e-15) {
  V <- t(V)
  W <- t(W)
  H <- t(H)
  W <- ccd_update_factors_rcpp(V,H,W,e)
  return(t(W))
}
