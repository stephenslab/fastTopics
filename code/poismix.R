# Maximize a Poisson likelihood in which the Poisson rate for the jth
# sample is r[j] = sum(L[j,]*x).
fitpoismix <- function (L, w, x, numiter = 100,
                        qp.solver = c("quadprog","activeset"),
                        beta = 0.75, suffdecr = 0.01, minstepsize = 1e-10,
                        e = 1e-8, delta = 1e-10, verbose = TRUE) {
  qp.solver <- match.arg(qp.solver)
    
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
    out <- fitpoismix.update(L,w,x,f,qp.solver,e,delta,beta,suffdecr,
                             minstepsize)
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
fitpoismix.update <- function (L, w, x, f,
                               qp.solver = c("quadprog","activeset"),
                               e = 1e-8, delta = 1e-6, beta = 0.75,
                               suffdecr = 0.01, minstepsize = 1e-10) {
  m         <- length(x)
  qp.solver <- match.arg(qp.solver)

  # Compute the gradient (g) and Hessian (H) at the current iterate.
  u <- drop(L %*% x) + e
  g <- drop((1 - w/u) %*% L)
  H <- as.matrix(crossprod((sqrt(w)/u)*L)) + delta*diag(m)
  
  # Compute a search direction, p, by minimizing p'*H*p/2 + p'*g,
  # where g is the gradient and H is the Hessian, subject to all
  # elements of x + p being non-negative. Although rather than solve
  # this problem directly, we instead minimize y'*H*y/2 + y'*(g - H*x)
  # subject to y being non-negative, then set p = y - x.
  ghat <- drop(g - H %*% x)
  if (qp.solver == "quadprog") {
    out <- quadprog::solve.QP(H,-ghat,Matrix::Diagonal(m))
    p   <- out$solution - x
  } else if (qp.solver == "activeset") {
    z <- activeset(H,ghat,x)
    p <- z - x
  }

  # Perform backtracking line search to determine a suitable step
  # size.
  a <- 0.99
  while (TRUE) {
    y <- x + a*p
    if (all(y >= 0)) {
      fnew <- cost.poismix(L,w,y,e)
      if (fnew <= f + suffdecr*a*dot(p,g))
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
