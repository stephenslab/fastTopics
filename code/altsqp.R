# Compute maximum-likelihood estimates for the Poisson topic model;
# equivalently, find a non-negative matrix factorization X = L*F' that
# optimizes the beta divergence objective.
altsqp <- function (X, F, L, numiter = 100, control = list(), verbose = TRUE) {

  # Get the optimization settings.
  control <- modifyList(altsqp_control_defaults,control,keep.null = TRUE)
  nc      <- control$nc
  e       <- control$e
  
  # Get the number of rows (n) and columns (m) of the counts matrix.
  n <- nrow(X)
  m <- ncol(X)

  # Set up the data structure to record the algorithm's progress.
  progress <- data.frame(iter      = 1:numiter,
                         objective = 0,
                         max.diff  = 0,
                         timing    = 0)
  
  # Iteratively apply the EM And SQP updates.
  if (verbose)
    cat("iter         objective max.diff\n")
  for (iter in 1:numiter) {

    # Save the current estimates of the factors and loadings.
    F0 <- F
    L0 <- L

    timing <- system.time({

      # Update the loadings ("activations").
      if (nc == 1)
        L <- altsqp.update.loadings(X,F,L,control)
      else {
        rows <- splitIndices(n,nc)
        L <- mclapply(rows,function (i) altsqp.update.loadings(X[i,],F,L[i,],
                                                               control),
                      mc.set.seed = FALSE,mc.allow.recursive = FALSE,
                      mc.cores = nc)
        L <- do.call(rbind,L)
        L[unlist(rows),] <- L
      }
          
      # Update the factors ("basis vectors").
      if (nc == 1)
        F <- altsqp.update.factors(X,F,L,control)
      else {
        cols <- splitIndices(m,nc)
        F <- mclapply(cols,function (j) altsqp.update.factors(X[,j],F[j,],L,
                                                              control),
                      mc.set.seed = FALSE,mc.allow.recursive = FALSE,
                      mc.cores = nc)
        F <- do.call(rbind,F)
        F[unlist(cols),] <- F
      }
    })

    # Compute the value of the objective (cost) function at the
    # current estimates of the factors and loadings.
    f <- cost(X,tcrossprod(L,F),e)
    d <- max(abs(tcrossprod(L,F) - tcrossprod(L0,F0)))
    progress[iter,"objective"] <- f
    progress[iter,"max.diff"]  <- d
    progress[iter,"timing"]    <- timing["elapsed"]
    if (verbose)
      cat(sprintf("%4d %+0.10e %0.2e\n",iter,f,d))
  }
  
  return(list(F = F,L = L,value = f,progress = progress))
}

altem <- function (X, F, L, numiter = 100, control = list()) {
  control <- modifyList(altsqp_control_defaults,control,keep.null = TRUE)
    
  # Get the number of rows (n) and columns (m) of the counts matrix.
  n <- nrow(X)
  m <- ncol(X)

 # Run DAAREM.
 out <- suppressWarnings(
   fpiter(log(c(as.vector(F),as.vector(L))),daarem.update,daarem.objective,
          control = list(maxiter = numiter,tol = 0),
          X = X,control.update = control))
  vars  <- project.iterate(out$par,n,m)
  value <- out$objfn.track[-1]
  return(list(F = vars$F,L = vars$L,value = value))
}

# These are the default optimization settings used in altsqp.
altsqp_control_defaults <- c(list(nc = 1),mixsqp_control_defaults)

# Update all the loadings with the factors remaining fixed.
altsqp.update.loadings <- function (X, F, L, control) {
  n <- nrow(X)
  for (i in 1:n) {
    L[i,] <- altsqp.update.em(F,X[i,],L[i,],control$e)
    L[i,] <- altsqp.update.sqp(F,X[i,],L[i,],control)
  }
  return(L)  
}

# Update all the factors with the loadings remaining fixed.
altsqp.update.factors <- function (X, F, L, control) {
  m <- ncol(X)
  for (j in 1:m) {
    F[j,] <- altsqp.update.em(L,X[,j],F[j,],control$e)
    F[j,] <- altsqp.update.sqp(L,X[,j],F[j,],control)
  }
  return(F)
}

# Run one EM update for the alternating SQP method.
altsqp.update.em <- function (B, w, y, e) {

  # Remove any counts that are exactly zero.
  ws <- sum(w)
  bs <- colSums(B)
  i  <- which(w > 0)
  w  <- w[i]
  B  <- B[i,]

  # Run an EM update for the modified problem.
  out <- mixem(scale.cols(B,ws/bs),w/ws,y*bs/ws,1,e)

  # Recover the solution to the unmodified problem.
  return(out$x*ws/bs)
}

# Run one SQP update for the alternating SQP method.
altsqp.update.sqp <- function (B, w, y, control) {
    
  # Remove any counts that are exactly zero.
  ws <- sum(w)
  bs <- colSums(B)  
  i  <- which(w > 0)
  w  <- w[i]
  B  <- B[i,]

  # Run an SQP update for the modified problem.
  out <- mixsqp(scale.cols(B,ws/bs),w/ws,y*bs/ws,1,control)

  # Recover the solution to the unmodified problem.
  return(out$x*ws/bs)
}

# Retrieve the m x k matrix of factors (F) and n x k matrix of
# loadings (L), and project them onto the non-negative quadrant.
project.iterate <-  function (vars, n, m) {
  nv <- length(vars)
  k  <- nv/(n + m)  
  F  <- matrix(vars[seq(1,m*k)],m,k)
  L  <- matrix(vars[seq(m*k+1,nv)],n,k)
  return(list(F = exp(F),L = exp(L)))
}

# This implements the objfn argument for the fpiter and daarem calls
# above.
daarem.objective <- function (vars, X, control.update) {
  n    <- nrow(X)
  m    <- ncol(X)
  vars <- project.iterate(vars,n,m)
  return(cost(X,tcrossprod(vars$L,vars$F),control.update$e))
}

# This implements the fixptfn argument for the fpiter and daarem
# calls above.
daarem.update <- function (vars, X, control.update) {
  n    <- nrow(X)
  m    <- ncol(X)
  vars <- project.iterate(vars,n,m)
  F    <- vars$F
  L    <- vars$L

  # Update the loadings ("activations").
  L <- altsqp.update.loadings(X,F,L,control.update)
  
  # Update the factors ("basis vectors").
  F <- altsqp.update.factors(X,F,L,control.update)

  cat(sprintf("%+0.10e\n",cost(X,tcrossprod(L,F),control.update$e)))
  return(log(c(as.vector(F),as.vector(L))))
}
