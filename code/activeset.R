# Active-set method for solving the convex quadratic program
#
#   minimize x'*H*x/2 + y'*g
#   subject to x >= 0.
#
# This implementation closely follows Algorithm 16.3 from Nocedal &
# Wright (2006).
activeset <- function (H, g, x0, maxiter = min(100,length(g) + 1),
                       convtol = 1e-15, zerosearchdir = 1e-15) {

  # Initialize the solution, and get the working set.
  n <- length(x0)
  x <- x0
  S <- (x >= 0)
  
  # Repeat until we have reached the maximum number of iterations, or
  # until the convergence criteria have been met.
  for (iter in 1:maxiter) {
      
    # Compute the search direction by solving the equality-constrained
    # subproblem.
    b    <- drop(H %*% x + g)
    i    <- which(S)
    p    <- rep(0,n)
    if (length(i) > 0)
      p[i] <- tryCatch(drop(solve(H[i,i],-b[i])),error = function (e) {
        browser()
      })

    # Check that the search direction is close to zero.
    if (max(abs(p)) < zerosearchdir) {
      j <- which(!S)

      # If the Lagrange multipliers for all co-ordinates in the
      # working set are positive, we have reached the solution.
      if (length(j) == 0 | all(b[j] >= -convtol))
        break

      # Find the co-ordinate with the smallest Lagrange multiplier,
      # and remove it from the working set.
      S[j[which.min(b[j])]] <- TRUE
    } else {

      # Determine the step size, check for any blocking constraints,
      # and adjust the working set if there are any blocking
      # constraints.
      a <- 1
      i <- which(S & p < 0)
      if (length(i) > 0) {

        # Check for a blocking constraint.
        y <- -x[i]/p[i]
        j <- which.min(y)
        if (y[j] < 1) {

          # There is a blocking constraint; adjust the step size to
          # retain feasibility, and add this blocking constraint to
          # the working set.
          a       <- y[j]
          S[i[j]] <- FALSE
        }
      }

      # Move to the new iterate along the search direction.
      x <- x + a*p
    }
  }

  return(x)
}
