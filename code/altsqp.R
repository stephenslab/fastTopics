# TO DO: Explain here what this function does, and how to use it.
#
# TO DO:
#
#   - Accommodate case when some weights (w) are zero.
#
#   - Implement function "fitpoismix.update" that implements a single
#     iteration of fitpoismix.
#
fitpoismix <- function (L, w, x, numiter = 1000, beta = 0.75, suffdecr = 0.01,
                        minstepsize = 1e-10, e = 1e-8, verbose = TRUE) {

  # Initialize the output data frame.
  progress <- data.frame(iter      = 1:numiter,
                         objective = 0,
                         max.diff  = 0,
                         step.size = 0)

  # Compute the value of the objective at the initial estimate.
  f <- cost.poismix(L,w,x,e)
  
  # Repeat until we reach the number of requested iterations.
  if (verbose)
    cat("iter         objective max.diff step.size\n")
  for (i in 1:numiter) {

    # Save the current estimates of the Poisson rates.
    x0 <- x

    # Compute the gradient (g) and Hessian (H) at the current iterate.
    u <- 1/(drop(L %*% x) + e)
    g <- drop((1 - w/u) %*% L)
    H <- crossprod((sqrt(w)*u)*L)
    
    # Compute a search direction, p, by minimizing p'*H*p/2 + p'*g,
    # where g is the gradient and H is the Hessian, subject to all
    # elements of x + p being non-negative.
    # TO DO.

    # Perform backtracking line search to determine a suitable step
    # size.
    a <- 0.99
    while (TRUE) {
      y    <- x + a*p
      fnew <- cost.poismix(L,w,y,e)
      if (fnew <= f + a*suffdecr*dot(p,g))
        break
      a <- a * beta
      if (a < minstepsize)
        break
    }

    # Move to the new iterate.
    x <- y
    f <- fnew
  }

  return(list(x = x,value = f,progress = progress))
}

# Return the cost function minimized by fitpoismix; it is the negative
# Poisson likelihood in which the Poisson rate for each sample is
# given by the matrix-vector product L*x.
cost.poismix <- function (L, w, x, e = 1e-8) {
 y <- drop(L %*% x) + e
 if (all(y > 0))
   return(sum(y - w*log(y)))
 else
   return(Inf)
}
