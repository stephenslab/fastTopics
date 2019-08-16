# Compute maximum-likelihood estimates of the mixture proportions in a
# mixture model by iterating the SQP updates for a fixed number of
# 
#' @importFrom utils modifyList
#' @importFrom Rcpp evalCpp
#' 
mixsqp <- function (L, w, x0, numiter, control = list(), verbose = FALSE) {

  # Get the number of rows (n) and columns (m) of the matrix L.
  n <- nrow(L)
  m <- ncol(L)

  # Get the optimization settings.
  control <- modifyList(mixsqp_control_default(),control,keep.null = TRUE)
  e       <- control$e
  
  # Run the updates implemented in C++.
  if (verbose)
    cat("iter objective function\n")
  control$maxiteractiveset <- m + 1
  x <- drop(mixsqp_rcpp(L,w,x0,rep(e,n),numiter,control,verbose))
  
  # Return (1) the estimate of the solution and (2) the value of the
  # objective at this estimate.
  return(list(x = x,value = mixobjective(L,w,x,e)))
}

# These are the default optimization settings used in the "mixsqp"
# function.
mixsqp_control_default <- function()
  list(activesetconvtol = 1e-10,
       zerothreshold    = 1e-10,
       zerosearchdir    = 1e-15,
       suffdecr         = 0.01,
       stepsizereduce   = 0.75,
       minstepsize      = 1e-10,
       e                = 1e-15)

# Compute the value of the mix-SQP objective at x; arguments L and w
# specify the objective, and e is a vector in which the entries can be
# set to small, positive numbers, or to zero.
mixobjective <- function (L, w, x, e) {
 y <- drop(L %*% x) + e
 if (all(y > 0))
   return(sum(x) - sum(w * log(y)))
 else
   return(Inf)
}
