# Maximize a Poisson likelihood in which the Poisson rate for the jth
# sample is r[j] = sum(L[j,]*x).
fitpoismix <- function (L, w, x, numiter = 100, beta = 0.75, suffdecr = 0.01,
                        minstepsize = 1e-10, e = 1e-8, verbose = TRUE) {

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
    out <- fitpoismix.update(L,w,x,f,e,suffdecr,minstepsize)
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

# Implement a single iteration of fixpoismix.
fitpoismix.update <- function (L, w, x, f, e, suffdecr, minstepsize) {
  m <- length(x)
    
  # Compute the gradient (g) and Hessian (H) at the current iterate.
  u <- drop(L %*% x) + e
  g <- drop((1 - w/u) %*% L)
  H <- crossprod((sqrt(w)/u)*L)
    
  # Compute a search direction, p, by minimizing p'*H*p/2 + p'*g,
  # where g is the gradient and H is the Hessian, subject to all
  # elements of x + p being non-negative.
  out <- solve.QP(H,-g,Diagonal(m),-x)
  p   <- out$solution
  
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

  # Output the new iterate (x), the value of the objective at this
  # point (f), and the step size (a).
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

