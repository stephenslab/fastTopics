# Compute maximum-likelihood estimates of the mixture proportions in a
# mixture model by iterating the SQP updates for a fixed number of
# iterations. This is mainly used for testing the C++ implementation.
# See the comments attached to the "mixem" C++ function for an
# explanation of the inputs.
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
  control$maxiteractiveset <- m + 1
  x <- drop(mixsqp_rcpp(L,w,x0,rep(e,n),numiter,control,verbose))
  
  # Return the updated estimate of the mixture weights.
  return(x)
}

# Returns the default algorithm settings for mixsqp.
mixsqp_control_default <- function()
  list(activesetconvtol = 1e-10,
       zerothreshold    = 1e-10,
       zerosearchdir    = 1e-15,
       suffdecr         = 0.01,
       stepsizereduce   = 0.75,
       minstepsize      = 1e-10,
       e                = 1e-15)
