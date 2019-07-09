#' @title Alternating SQP method for non-negative matrix factorization
#'   and topic modeling.
#' 
#' @description Compute maximum-likelihood estimates for the Poisson
#'   topic model; equivalently, find a non-negative matrix factorization
#'   X = L*F' that optimizes the beta divergence objective.
#'
#' @param X Describe X here.
#'
#' @param F Describe F here.
#'
#' @param L Describe L here.
#'
#' @param numiter Describe numiter here.
#'
#' @param control Describe control here.
#'
#' @param verbose Describe verbose here.
#'
#' @importFrom utils modifyList
#' @importFrom parallel splitIndices
#' @importFrom parallel mclapply
#' 
#' @export
#' 
altsqp <- function (X, F, L, numiter = 100, control = list(), verbose = TRUE) {

  # Get the optimization settings.
  control <- modifyList(altsqp_control_default(),control,keep.null = TRUE)
  b0      <- control$b0
  exiter0 <- control$exiter0
  bmaxinc <- control$bmaxinc
  binc    <- control$binc
  bred    <- control$bred
  nc      <- control$nc
  e       <- control$e
  
  # Get the number of rows (n) and columns (m) of the counts matrix.
  n <- nrow(X)
  m <- ncol(X)

  # Compute the value of the objective at the initial iterate.
  f     <- cost(X,tcrossprod(L,F),e)
  fbest <- f
  
  # These are additional quantities used to implement the
  # extrapolation scheme of Ang & Gillis (2019).
  bmax  <- 0.99
  b     <- 0
  Fy    <- F
  Fbest <- F
  Ly    <- L 
  Lbest <- L
  
  # Set up the data structure to record the algorithm's progress.
  progress <- data.frame(iter      = 1:numiter,
                         objective = 0,
                         max.diff  = 0,
                         beta      = 0,
                         timing    = 0)

  # Iteratively apply the EM And SQP updates.
  if (verbose)
    cat("iter         objective max.diff    beta\n")
  for (iter in 1:numiter) {

    # Store the value of the objective at the current iterate.
    f0 <- f

    # When the time is right, initiate the extrapolation scheme.
    if (iter == exiter0)
      b <- b0
      
    timing <- system.time({

      # Update the factors ("basis vectors").
      if (nc == 1)
        Fn <- altsqp.update.factors(X,Fy,Ly,control)
      else {
        cols <- splitIndices(m,nc)
        Fn <- mclapply(cols,
                function (j) altsqp.update.factors(X[,j],Fy[j,],Ly,control),
                mc.set.seed = FALSE,mc.allow.recursive = FALSE,mc.cores = nc)
        Fn <- do.call(rbind,Fn)
        Fn[unlist(cols),] <- Fn
      }

      # Compute the extrapolated update for the factors.
      Fy <- pmax(Fn + b*(Fn - F),0)
      
      # Update the loadings ("activations").
      if (nc == 1)
        Ln <- altsqp.update.loadings(X,Fy,Ly,control)
      else {
        rows <- splitIndices(n,nc)
        Ln <- mclapply(rows,
               function (i) altsqp.update.loadings(X[i,],Fy,Ly[i,],control),
               mc.set.seed = FALSE,mc.allow.recursive = FALSE,mc.cores = nc)
        Ln <- do.call(rbind,Ln)
        Ln[unlist(rows),] <- Ln
      }

      # Compute the extrapolated update for the loadings.
      Ly <- pmax(Ln + b*(Ln - L),0)
    })

    # Compute the value of the objective (cost) function at the
    # extrapolated solution for the factors (F) and the basic co-ordinate
    # ascent solution for the loadings (L).
    f <- cost(X,tcrossprod(Ln,Fy),e)

    if (b == 0) {
      F <- Fn
      L <- Ln
    } else {

      # Before updating the extrapolation parameter, store the current
      # value.
      b0 <- b
      if (f > f0) {

        # The solution did not improve, so restart the extrapolation scheme.
        Fy   <- F
        Ly   <- L
        bmax <- b0
        b    <- bred*b
      } else {
        
        # The solution is improved; keep the basic co-ordinate ascent
        # update as well.
        F    <- Fn
        L    <- Ln
        b    <- min(bmax,binc * b)
        bmax <- min(0.99,bmaxinc * bmax)
      }
    }        

    # If the solution is improved, update the current best solution.
    if (f < f0) {
      d     <- max(abs(tcrossprod(Lbest,Fbest) - tcrossprod(Ln,Fy)))
      Fbest <- Fy
      Lbest <- Ln 
      fbest <- f
    } else
      d <- 0

    # Record the algorithm's progress.
    #
    # TO DO: Add comments here.
    # 
    progress[iter,"objective"] <- fbest
    progress[iter,"max.diff"]  <- d
    progress[iter,"beta"]      <- b
    progress[iter,"timing"]    <- timing["elapsed"]
    if (verbose)
      cat(sprintf("%4d %+0.10e %0.2e %0.1e\n",iter,fbest,d,b))
  }

  # TO DO: Add comments here.
  return(list(F = Fbest,L = Lbest,value = fbest,progress = progress))
}

#' @rdname altsqp
#' 
#' @export
#' 
altsqp_control_default <- function()
  c(mixsqp_control_default(),
    list(nc = 1,exiter0 = 10,b0 = 0.5,bmaxinc = 1.05,binc = 1.1,bred = 0.75))

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
  i  <- which(w > 0)
  w  <- w[i]
  B  <- B[i,]
  ws <- sum(w)
  bs <- colSums(B)

  # Run an EM update for the modified problem.
  out <- mixem(scale.cols(B,ws/bs),w/ws,y*bs/ws,1,e)

  # Recover the solution to the unmodified problem.
  return(out$x*ws/bs)
}

# Run one SQP update for the alternating SQP method.
altsqp.update.sqp <- function (B, w, y, control) {
    
  # Remove any counts that are exactly zero.
  i  <- which(w > 0)
  w  <- w[i]
  B  <- B[i,]
  ws <- sum(w)
  bs <- colSums(B)  
  
  # Run an SQP update for the modified problem.
  out <- mixsqp(scale.cols(B,ws/bs),w/ws,y*bs/ws,1,control)

  # Recover the solution to the unmodified problem.
  return(out$x*ws/bs)
}
