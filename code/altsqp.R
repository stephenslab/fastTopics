# TO DO: Explain here what this function does, and how to use it.
altsqp <- function (X, F, L, numiter = 1000, e = 1e-8, verbose = TRUE) {
  n <- nrow(X)
  m <- ncol(X)
  
  progress <- data.frame(iter = 1:numiter,objective = 0,
                         max.diff = 0,timing = 0)

  # Repeat until we reach the number of requested iterations.
  if (verbose)
    cat("iter         objective max.diff\n")
  for (iter in 1:numiter) {

    # Save the current estimates of the factors and loadings.
    F0 <- F
    L0 <- L

    timing <- system.time({
        
      # Update the loadings ("activations").
      for (i in 1:n) {
        fi    <- cost.poismix(F,X[i,],L[i,],e)
        out   <- fitpoismix.update(F,X[i,],L[i,],fi,e)
        L[i,] <- out$x
      }

      # Update the factors ("basis vectors").
      for (j in 1:m) {
        fj    <- cost.poismix(L,X[,j],F[j,],e)
        out   <- fitpoismix.update(L,X[,j],F[j,],fj,e)
        F[j,] <- out$x
      }
    })

    # Compute the value of the objective (cost) function at the
    # current estimates of the factors and loadings.
    f <- cost(X,tcrossprod(L,F),e)
    d <- max(max(abs(F - F0)),max(abs(L - L0)))
    progress[iter,"objective"] <- f
    progress[iter,"max.diff"]  <- d
    progress[iter,"timing"]    <- timing["elapsed"]
    if (verbose)
      cat(sprintf("%4d %+0.10e %0.2e\n",iter,f,d))
  }

  return(list(F = F,L = L,value = f,progress = progress))
}

# Maximize a Poisson likelihood in which the Poisson rate for the jth
# sample is r[j] = sum(L[j,]*x).
fitpoismix <- function (L, w, x, numiter = 100, beta = 0.75, suffdecr = 0.01,
                        minstepsize = 1e-10, e = 1e-8, delta = 1e-10,
                        verbose = TRUE) {

  # Remove rows with zero weights.
  rows <- which(w > 0)
  w    <- w[rows]
  L    <- L[rows,]

  # Initialize the output data frame.
  progress <- data.frame(iter      = 1:numiter,
                         objective = 0,
                         max.diff  = 0,
                         step.size = 0)

  # Compute the value of the objective at the initial solution.
  f <- cost.poismix(L,w,x,e)
  
  # Repeat until we reach the number of requested iterations.
  if (verbose)
    cat("iter           objective max.diff step.size\n")
  for (i in 1:numiter) {
    x0  <- x
    out <- fitpoismix.update(L,w,x,f,e,delta,beta,suffdecr,minstepsize)
    x   <- out$x
    f   <- out$f
    a   <- out$a
    d   <- max(abs(x - x0))
    progress[i,"objective"] <- f
    progress[i,"max.diff"]  <- d
    progress[i,"step.size"] <- a
    if (verbose)
      cat(sprintf("%4d %+0.12e %0.2e %0.3e\n",i,f,d,a))
  }

  return(list(x = x,value = f,progress = progress))
}

# Implements a single iteration of fixpoismix.
fitpoismix.update <- function (L, w, x, f, e = 1e-8, delta = 1e-6, beta = 0.75,
                               suffdecr = 0.01, minstepsize = 1e-10) {
  m <- length(x)
    
  # Compute the gradient (g) and Hessian (H) at the current iterate.
  u <- drop(L %*% x) + e
  g <- drop((1 - w/u) %*% L)
  H <- as.matrix(crossprod((sqrt(w)/u)*L)) + delta*diag(m)
  
  # Compute a search direction, p, by minimizing p'*H*p/2 + p'*g,
  # where g is the gradient and H is the Hessian, subject to all
  # elements of x + p being non-negative.
  out <- solve.QP(H,-g,Diagonal(m),-x)
  p   <- out$solution

  # Perform backtracking line search to determine a suitable step
  # size.
  a <- 0.99
  while (TRUE) {
    y <- x + a*p
    if (all(y >= 0)) {
      fnew <- cost.poismix(L,w,y,e)
      if (fnew <= f + a*suffdecr*dot(p,g))
        break
    }
    if (a*beta < minstepsize)
      break
    a <- a * beta
  }

  # Output the new iterate (x), the value of the objective at this
  # point (f), and the step size (a).
  fnew <- cost.poismix(L,w,y,e)
  return(list(x = y,f = fnew,a = a))
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
