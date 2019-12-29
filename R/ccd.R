# TO DO: Explain here what this function does, and how to use it.
ccd_update_factors <- function (X, A, B, e = 1e-15) {
  B1 <- B + 0
  AB <- A %*% B
  A  <- t(A)
  ccd_update_factors_rcpp(X,A,B1,AB,e)
  return(B1)
}

# TO DO: Explain here what this function does, and how to use it.
ccd_update_loadings <- function (X, A, B, e = 1e-15) {
  X  <- t(X)
  AB <- t(A %*% B)
  A  <- t(A)
  # ccd_update_loadings_rcpp(X,A,B,AB,e)
  ccd_update_factors_rcpp(X,B,A,AB,e)
  return(t(A))
}

ccd <- function (X, A, B, e = 1e-15) {
  m  <- ncol(X)
  A  <- t(A)
  ab <- rep(0,m)
  v  <- rep(0,m)
  AB <- t(A) %*% B
  ccd_rcpp(X,A,B,AB,ab,v,e)
  A  <- t(A)
  return(list(A = A,B = B))
}
