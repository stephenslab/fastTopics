# Compute maximum-likelihood estimates for the Poisson topic model;
# equivalently, find a non-negative matrix factorization X = L*F' that
# optimizes the beta divergence objective.
altsqp <- function (X, F, L, numiter = 100, method = c("fpiter","daarem"),
                    control = list(), verbose = TRUE) {

  # Get the optimization settings.
  method  <- match.arg(method)
  control <- modifyList(altsqp_control_defaults,control,keep.null = TRUE)

  # Get the number of rows (n) and columns (m) of the counts matrix.
  n <- nrow(X)
  m <- ncol(X)

  # Iteratively apply the EM And SQP updates.
  if (verbose)
      cat("         objective max.diff time(s)\n")
  if (method == "fpiter")
    out <- suppressWarnings(
      fpiter(c(as.vector(F),as.vector(L)),altsqp.update,altsqp.objective,
             list(maxiter = numiter,order = control$order,tol = 0,
                  mon.tol = 0.05,kappa = 20,alpha = 1.2),
             X = X,control.update = control,verbose = verbose))
  else {
    # TO DO.
  }

  # Return (1) the estimated factors, (2) the estimated loadings, and
  # (3) the value of the objective at each iteration.
  vars <- project.iterate(out$par)
  return(list(F = vars$F,L = vars$L,value = out$objfn.track[-1]))
}

# These are the default optimization settings used in altsqp.
altsqp_control_defaults <- c(list(nem = 1,nsqp = 4,nc = 1,order = 10),
                             mixsqp_control_defaults)
       
# Retrieve the m x k matrix of factors (F) and n x k matrix of
# loadings (L), and project them onto the non-negative quadrant.
project.iterate <-  function (vars, n, m) {
  nv <- length(vars)
  k  <- nv/(n + m)  
  F  <- matrix(vars[seq(1,m*k)],p,k)
  L  <- matrix(vars[seq(m*k+1,nv)],k)
  return(list(F = pmax(F,0),L = pmax(L,0)))
}

# This implements the objfn argument for the fpiter and daarem calls
# above.
altsqp.objective <- function (vars, X, control.update, verbose) {
  n    <- nrow(X)
  m    <- ncol(X)
  e    <- control.update$e
  vars <- project.iterate(vars,n,m)
  return(cost(X,tcrossprod(vars$L,vars$F),e))
}

# This implements the fixptfn argument for the fpiter and daarem
# calls above.
altsqp.update <- function (vars, X, control.update, verbose) {
  n    <- nrow(X)
  m    <- ncol(X)
  e    <- control.update$e
  vars <- project.iterate(vars,n,m)
  F    <- vars$F
  L    <- vars$L
  
  timing <- system.time({

    # Update the loadings ("activations").
    if (nc == 1)
      L <- altsqp.update.loadings(X,F,L,control.update)
    else {
      rows <- splitIndices(n,nc)
      L <- mclapply(rows,function (i) altsqp.update.loadings(X[i,],F,L[i,],
                                                             control.update))
      L <- do.call(rbind,L)
      L[unlist(rows),] <- L
    }
          
    # Update the factors ("basis vectors").
    if (nc == 1)
      F <- altsqp.update.factors(X,F,L,control.update)
    else {
      cols <- splitIndices(m,nc)
      F <- mclapply(cols,function (j) altsqp.update.factors(X[,j],F[j,],L,
                                                            control.update))
      F <- do.call(rbind,F)
      F[unlist(cols),] <- F
    }
  })

  # Report the algorithm's progress, if requested.
  if (verbose) {
    f <- cost(X,tcrossprod(L,F),e)
    d <- max(max(abs(F/rowMeans(F) - F0/rowMeans(F0))),
             max(abs(L/rowMeans(L) - L0/rowMeans(L0))))
    cat(sprintf("%+0.10e %0.2e %0.2f\n",f,d,timing["elapsed"]))
  }
}

# Update all the loadings with the factors remaining fixed.
altsqp.update.loadings <- function (X, F, L, control) {
  n <- nrow(X)
  for (i in 1:n) {
    if (nem > 0)
      L[i,] <- altsqp.update.em(F,X[i,],L[i,],control$nem,control$e)
    if (nsqp > 0)
      L[i,] <- altsqp.update.sqp(F,X[i,],L[i,],control$nsqp,control)
  }
  return(L)  
}

# Update all the factors with the loadings remaining fixed.
altsqp.update.factors <- function (X, F, L, control) {
  m <- ncol(X)
  for (j in 1:m) {
    if (nem > 0)
      F[j,] <- altsqp.update.em(L,X[,j],F[j,],control$nem,control$e)
    if (nsqp > 0)
      F[j,] <- altsqp.update.sqp(L,X[,j],F[j,],control$nsqp,control)
  }
  return(F)
}

# Run a fixed number of EM updates for the alternating SQP method.
altsqp.update.em <- function (B, w, y, numiter, e) {

  # Remove any counts that are exactly zero.
  ws <- sum(w)
  bs <- colSums(B)
  i  <- which(w > 0)
  w  <- w[i]
  B  <- B[i,]

  # Run the EM updates for the modified problem.
  out <- mixem(scale.cols(B,ws/bs),w/ws,y*bs/ws,numiter,e)

  # Recover the solution to the unmodified problem.
  return(out$x*ws/bs)
}

# Run a fixed number of SQP updates for the alternating SQP method.
altsqp.update.sqp <- function (B, w, y, numiter, control) {
    
  # Remove any counts that are exactly zero.
  ws <- sum(w)
  bs <- colSums(B)  
  i  <- which(w > 0)
  w  <- w[i]
  B  <- B[i,]

  # Run the SQP updates for the modified problem.
  out <- mixsqp(scale.cols(B,ws/bs),w/ws,y*bs/ws,numiter,control)
  
  # Recover the solution to the unmodified problem.
  return(out$x*ws/bs)
}
