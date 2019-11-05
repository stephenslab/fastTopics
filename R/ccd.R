#' @keywords internal
#'
#' @export
#' 
ccd <- function (X, A, B, numiter = 100, e = 1e-15) {

  # This variable is used to keep track of the algorithm's progress;
  # it stores the value of the cost function at each iteration.
  value <- rep(0,numiter)

  # Execute the CCD updates.
  A  <- t(A)
  n  <- ncol(X)
  ab <- rep(0,m)
  v  <- rep(0,m)
  for (i in 1:numiter) {
    AB <- t(A) %*% B
    ccd_rcpp(X,A,B,AB,ab,v,e)
    value[i] <- cost(X,t(A),B,e)
  }
  A <- t(A)
  
  return(list(A = A,B = B,value = value))
}
