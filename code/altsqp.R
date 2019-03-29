# TO DO: Explain here what this function does, and how to use it.
altsqp <- function (X, F, L, numiter = 1000, e = 1e-8, nc = 1,
                    verbose = TRUE) {
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
      if (nc == 1)  
        L <- altsqp.update.loadings(X,F,L,e)
      else {
        rows <- splitIndices(n,nc)
        L <- mclapply(rows,function(i) altsqp.update.loadings(X[i,],F,L[i,],e))
        L <- do.call(rbind,L)
        L[unlist(rows),] <- L
      }

      # Update the factors ("basis vectors").
      if (nc == 1)
        F <- altsqp.update.factors(X,F,L,e)
      else {
        cols <- splitIndices(m,nc)
        F <- mclapply(cols,function(j) altsqp.update.factors(X[,j],F[j,],L,e))
        F <- do.call(rbind,F)
        F[unlist(cols),] <- F
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

# This implements the loadings update step in the alternating SQP method.
altsqp.update.loadings <- function (X, F, L, e) {
  n <- nrow(X)
  for (i in 1:n) {
    fi    <- cost.poismix(F,X[i,],L[i,],e)
    out   <- fitpoismix.update(F,X[i,],L[i,],fi,e = e)
    L[i,] <- out$x
  }
  return(L)
}

# This implements the factors update step in the alternating SQP method.
altsqp.update.factors <- function (X, F, L, e) {
  m <- ncol(X)
  for (j in 1:m) {
    fj    <- cost.poismix(L,X[,j],F[j,],e)
    out   <- fitpoismix.update(L,X[,j],F[j,],fj,e = e)
    F[j,] <- out$x
  }
  return(F)
}

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
  browser()
  if (qp.solver == "quadprog") {
    out <- quadprog::solve.QP(H,drop(H %*% x - g),Matrix::Diagonal(m))
    p   <- out$solution - x
  } else if (qp.solver == "activeset") {
    z    <- 1 * x
    ghat <- drop(g - H %*% x)
    activeset_rcpp(H,ghat,x,z,100,1e-10,0,1e-10)
    p <- z - x
  }

  # Perform backtracking line search to determine a suitable step
  # size.
  a <- 1
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
