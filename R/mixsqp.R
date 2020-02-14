# Compute maximum-likelihood estimates of the mixture proportions in a
# mixture model by iterating the mix-SQP updates for a fixed number of
# iterations. This is mainly used for testing the C++ implementation.
# See the comments attached to the "mixsqp" C++ function for an
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
  x <- drop(mixsqp_rcpp(L,w,x0,numiter,control)$x)
  
  # Return the updated estimate of the mixture weights.
  return(x)
}

# These are the default settings used for running mix-SQP.
mixsqp_control_default <- function()
  list(convtol.activeset         = 1e-10,
       zero.threshold.solution   = 1e-6,
       zero.threshold.searchdir  = 1e-10,
       suffdecr.linesearch       = 0.01,
       stepsizereduce            = 0.75,
       minstepsize               = 1e-8,
       identity.contrib.increase = 10,
       eps                       = 1e-8,
       maxiter.activeset         = 20)
