# TO DO: Explain here what this function does, and how to use it.
fitpoismix <- function (L, w, x, numiter = 1000, e = 1e-8) {

  # Repeat until we reach the number of requested iterations.
  for (i in 1:numiter) {

    # Save the current estimates of the Poisson rates.
    x0 <- x

    # Compute the gradient and Hessian at x.
    # TO DO.

    # Compute a search direction p by minimizing p'*H*p/2 + p'*g,
    # where g is the gradient and H is the Hessian, subject to all
    # elements of x + p being positive.
    # TO DO.

    # Perform backtracking line search to determine a suitable step
    # size.
    # TO DO.
  }
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
