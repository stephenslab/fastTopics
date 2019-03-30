# Active-set method for solving the convex quadratic program
#
#   minimize x'*H*x/2 + y'*g
#   subject to x >= 0.
#
activeset <- function (H, g, x0, maxiter = 100) {

  # Get the working set.
  n <- length(x0)
  x <- x0
  S <- (x >= 0)

  for (iter in 1:maxiter) {
      
    # Compute the search direction by solving the equality-constrained
    # subproblem.
    b    <- H %*% x + g
    i    <- which(S)
    p    <- rep(0,n)
    p[i] <- qr.solve(H[i,i],-b[i])

    # Check that the search direction is close to zero.
    if (max(abs(p)) < 1e-15) {
      i <- which(!S)

      # If the Lagrange multipliers for all co-ordinates in the
      # working set
      if (length(i) == 0)
        break
      else if (all(b[i] >= -1e-15))
        break
      else {

        # Find the co-ordinate with the smallest Lagrange multiplier,
        # and remove it from the working set.
        j    <- i[which.min(b[i])]
        S[j] <- TRUE
      }
    } else {

      # Determine the step size.
      a <- 1
      i <- which(S & p < 0)
      if (length(i) > 0) {
        y <- -x[i]/p[i]
        j <- which.min(y)

        # Check if there are any blocking constraints. If so, then the
        # step size will be less than 1.
        if (y[j] < 1) {
          a    <- y[j]
          j    <- i[j]
          S[j] <- FALSE
        }
      }

      # Move to the new iterate along the search direction.
      x <- x + a*p
    }
  }

  return(x)
}
