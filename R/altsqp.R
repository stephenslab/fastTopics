#' @title Alternating SQP Method for Optimizing Topic Models and Non-negative Matrix Factorizations
#' 
#' @description Compute maximum-likelihood estimates for the Poisson
#'   topic model with factors F and loadings L; equivalently, find a
#'   non-negative matrix factorization \eqn{L * F^T} of matrix X that
#'   optimizes the "generalized Kullback-Leibler" (or Bregman)
#'   divergence objective.
#'
#' @details Functions \code{loglik.poisson} and \code{loglik.multinom}
#' compute the log-likelihood for the Poisson and multinomial topic
#' models, excluding terms that do not depend on the model
#' parameters. These functions can be used to evaluate a solution
#' returned by \code{altsqp}.
#'
#' Use function \code{poisson2multinom} to obtain parameters for the
#' multinomial topic model from the Poisson topic model.
#' 
#' The \code{control} argument to \code{altsqp} is a list in which any
#' of the following named components will override the default
#' optimization algorithm settings (as they are defined by
#' \code{altsqp_control_default}):
#' 
#' \describe{
#' 
#' \item{\code{nc}}{This determines the \code{mc.cores} argument in
#'   calls to \code{\link[parallel]{mclapply}}.}
#'
#' \item{\code{numem}}{A non-negative number specifying the number of
#'   EM (\emph{i.e.}, multiplicative) updates to run at each outer loop
#'   iteration.}
#'
#' \item{\code{numsqp}}{A non-negative number specifying the number of
#'   SQP updates to run at each outer loop iteration.}
#' 
#' \item{\code{extrapolate}}{The iteration at which extrapolation is
#'   initiated. If \code{extrapolate > numiter}, extrapolation is not
#'   used.}
#'
#' \item{\code{beta.init}}{The initial setting of the extrapolation
#'   parameter.}
#'
#' \item{\code{beta.increase}}{The extrapolation parameter is
#'   increased by this amount whenever the update improves the
#'   solution. This is denoted by \eqn{\gamma} in Algorithm 3 of Ang &
#'   Gillis (2019).}
#'
#' \item{\code{beta.reduce}}{The extrapolation parameter is decreased
#'   by this amount whenever the update does not improve the solution.
#'   This is denoted by \eqn{\eta} in Algorithm 3 of Ang &
#'   Gillis (2019).}
#'
#' \item{\code{betamax.increase}}{The upper bound on the extrapolation
#'   parameter is increased by this amount whenever the update improves
#'   the solution. This is denoted by \eqn{\bar{\gamma}} in Algorithm 3
#'   of Ang & Gillis (2019).}
#' 
#' \item{\code{tol}}{A small, non-negative number specifying the
#'   convergence tolerance for the active-set step. Smaller values will
#'   result in higher quality search directions for the SQP algorithm
#'   but possibly a greater per-iteration computational cost.}
#'
#' \item{\code{zero.threshold}}{A small, non-negative number used to
#'   determine the "active set"; that is, it determines which entries of
#'   the solution are exactly zero. Any entries that are less than or
#'   equal to \code{zero.threshold} are considered to be exactly
#'   zero. Larger values of \code{zero.threshold} may lead to speedups
#'   for matrices with many columns, at the (slight) risk of prematurely
#'   zeroing some co-ordinates.}
#'
#' \item{\code{zero.searchdir}}{A small, non-negative number used to
#'   determine when the search direction in the active-set step is
#'   considered "small enough".}
#'
#' \item{\code{suffdecr}}{This parameter determines how stringent the
#'   "sufficient decrease" condition is for accepting a step size in the
#'   backtracking line search. Larger values will make the condition
#'   more stringent. This should be a positive number less than 1.}
#'
#' \item{\code{stepsizereduce}}{The multiplicative factor for
#'   decreasing the step size in the backtracking line search.}
#'
#' \item{\code{minstepsize}}{The smallest step size accepted by the
#'   line search step. Should be a number greater than 0 and at most 1.}
#'
#' \item{\code{e}}{The same as input argument \code{e} described above.}
#' }
#' 
#' @param X The n x m matrix of counts or pseudocounts. It can be a
#'   dense matrix or sparse matrix.
#'
#' @param fit A list containing two (dense, non-negative) matrices,
#'   \code{fit$F} and \code{fit$L}; for example, this can be the output
#'   of \code{altsqp}. The former is an m x k matrix of factors (or
#'   "basic vectors"), and the latter is an n x k matrix of loadings (or
#'   "activations"). For \code{loglik.multinom}, it is additionally
#'   required that each row of \code{fit$L} and each column of
#'    \code{fit$F} must sum to 1; \code{poisson2multinom} can be used to
#'   generate factors and loadings satisfying this requirement. For
#'   \code{altsqp}, this provides the initial estimates.
#'
#' @param numiter A positive integer specifying the number of
#'   updates to perform.
#'
#' @param version With \code{version = "R"}, a slower, but more easily
#'   debugged, implementation is used; with \code{version = "Rcpp"}, a
#'   faster, less easily debugged implementation is used.
#' 
#' @param control A list of parameters controlling the behaviour of
#'   the optimization algorithm. See \sQuote{Details}.
#'
#' @param verbose If \code{verbose = TRUE}, the algorithm's progress
#'   and a summary of the optimization settings are printed to the
#'   console.
#'
#' @param e A small, non-negative number that is added to the
#'   terms inside the logarithms to sidestep computing logarithms of
#'   zero. This prevents numerical problems at the cost of introducing a
#'   small (and typically very small) inaccuracy in the computation.
#'
#' @return \code{altsqp} returns a list object with the following
#' elements:
#'
#' \item{F}{A dense matrix containing estimates of the factors.}
#'
#' \item{L}{A dense matrix containing estimates of the loadings.}
#'
#' \item{value}{The value of the objective (or cost function) at the
#'   outputted values of \code{F} and \code{L}. The objective function
#'   is equal to \code{-loglik.poisson(X,out$,out$L)}, where \code{out}
#'   is the \code{altsqp} return value.}
#'
#' \item{progress}{A data frame containing more detailed information
#'   about the algorithm's progress. The data frame should have
#'   \code{numiter} rows. The columns of the data frame are: "iter", the
#'   SQP iteration; "objective", the value of the objective at the
#'   current estimate of the solution; "max.diff", the maximum
#'   difference in the Poisson rate parameters \code{tcrossprod(L,F)}
#'   between two the successive iterations; "beta", the current setting
#'   of the extrapolation parameter (zero means no extrapolation is
#'   used); and "timing", the elapsed time in seconds (based on
#'   \code{\link{system.time}}).}
#'
#' \code{poisson2multinom} returns an updated \code{fit}, in which the
#' factors \code{fit$F} and loadings \code{fit$L} are converted to
#' matrices of the same dimension containing, respectively, word
#' probabilities and topic probabilities.
#' 
#' \item{F}{A dense matrix containing estimates of the factors.}
#'
#' \item{L}{A dense matrix containing estimates of the loadings.}
#' 
#' \code{loglik.poisson} and \code{loglik.multinom} return a single
#' numeric value giving the value of the log-likelihood.
#' 
#' @references
#'
#' A. Ang and N. Gillis (2019). Accelerating nonnegative matrix
#' factorization algorithms using extrapolation. \emph{Neural Computation}
#' \bold{31}, 417-–439. \url{https://doi.org/10.1162/neco_a_01157}
#' 
#' @examples
#'
#' library(Matrix)
#' library(NNLM)
#'
#' # Generate a 300 x 400 data matrix to factorize. Less than 10% of the
#' # matrix elements should be nonzero.
#' set.seed(1)
#' n <- 300
#' m <- 400
#' k <- 3
#' F <- matrix(runif(m*k)/3,m,k)
#' L <- matrix(runif(n*k)/3,n,k)
#' X <- matrix(rpois(n*m,tcrossprod(L,F)),n,m)
#' X <- as(X,"dgCMatrix")
#' nnzero(X)/(n*m)
#' 
#' # Generate random initial estimates of the factors and loadings.
#' fit0 <- list(F = matrix(runif(m*k),m,k),
#'              L = matrix(runif(n*k),n,k))
#'
#' # Run 60 iterations of the sequential coordinate-wise descent
#' # algorithm implemented in the NNLM package. Note that nnmf does not
#' # accept a sparse matrix as input, so we need to provide it with a
#' # dense matrix instead.
#' fit1 <- suppressWarnings(
#'   nnmf(as.matrix(X),k,init = list(W = fit0$L,H = t(fit0$F)),
#'        method = "scd",loss = "mkl",max.iter = 60,rel.tol = 0, 
#'        inner.max.iter = 4,trace = 1,verbose = 0))
#'
#' # Run 60 coordinate-wise updates of the SQP method implemented in the
#' # fastTopics package.
#' fit2 <- altsqp(X,fit0,numiter = 60,version = "R",verbose = FALSE)
#' 
#' # Compare the Poisson log-likelihood at the two solutions; the
#' # likelihood should be higher at the the altsqp solution.
#' fit1$F <- t(fit1$H)
#' fit1$L <- fit1$W
#' print(loglik.poisson(X,fit1),digits = 14)
#' print(loglik.poisson(X,fit2),digits = 14)
#' 
#' # Compare the multinomial log-likelihood at the two solutions; again,
#' # the likelihood should be higher at the altsqp solution.
#' print(loglik.multinom(X,poisson2multinom(fit1)),digits = 14)
#' print(loglik.multinom(X,poisson2multinom(fit2)),digits = 14)
#' 
#' # Plot the improvement in the solution over time; the altsqp iterates
#' # (the solid, orange line) gets much closer to the best solution.
#' fbest    <- 31041.93896745
#' fit1$mkl <- n*m*fit1$mkl + sum(X - X*log(X + 1e-16))
#' plot(fit2$progress$iter,fit2$progress$objective - fbest,
#'      log = "y",col = "darkorange",type = "l",lwd = 2,xlab = "iteration",
#'      ylab = "distance from solution")
#' lines(1:60,fit1$mkl - fbest,col = "darkblue",lwd = 2,lty = "dashed")
#'
#' @useDynLib fastTopics
#' 
#' @importFrom Rcpp evalCpp
#' @importFrom utils modifyList
#' @importFrom Matrix t
#' @importFrom Matrix rowSums
#' @importFrom Matrix colSums
#' @importFrom Matrix mean
#' @importFrom parallel splitIndices
#' @importFrom parallel mclapply
#' 
#' @export
#' 
altsqp <- function (X, fit, numiter = 100, version = c("Rcpp", "R"),
                    control = list(), verbose = TRUE) {

  # Verify and process input matrix X. Each row and each column of the
  # matrix should have at least two positive entries.
  verify.matrix(X)
  if (!(all(rowSums(X > 0) >= 2) &
        all(colSums(X > 0) >= 2)))
    stop(paste("Each row and column of \"X\" should have at least two",
               "positive entries"))
  if (is.matrix(X) & mean(X > 0) < 0.1)
    message(paste("Input matrix \"X\" has less than 10% nonzero entries;",
                  "consider converting \"X\" to a sparse matrix to reduce",
                  "computational effort"))
  if (is.matrix(X) & is.integer(X))
    storage.mode(X) <- "double"

  # Input argument "fit" should be a list with elements "F" and "L".
  verify.fit(fit)
  F <- fit$F
  L <- fit$L
  
  # Verify and process input matrix F.
  verify.matrix(fit$F)
  if (!is.matrix(F))
    F <- as.matrix(F)
  if (is.integer(F))
    storage.mode(F) <- "double"
  if (all(F <= 0))
    stop("Some entries of input matrix \"fit$F\" should be positive")
  
  # Verify and process input matrix L.
  verify.matrix(fit$L)
  if (!is.matrix(L))
    L <- as.matrix(L)
  if (is.integer(L))
    storage.mode(L) <- "double"
  if (all(L <= 0))
    stop("Some entries of input matrix \"fit$L\" should be positive")

  # Check that matrices X, F and L are compatible.
  if (!(nrow(L) == nrow(X) &
        nrow(F) == ncol(X) &
        ncol(L) == ncol(F)))
    stop(paste("Dimensions of input arguments \"X\", \"fit$F\" and/or",
               "\"fit$L\ do not agree"))

  # Get the number of rows (n) and columns (m) of X, and the rank of
  # the matrix factorization (k).
  n <- nrow(X)
  m <- ncol(X)
  k <- ncol(F)

  # Determine which implementation to use.
  version <- match.arg(version)
  
  # Get the optimization settings.
  control <- modifyList(altsqp_control_default(),control,keep.null = TRUE)
  control$maxiteractiveset <- k + 1
  e       <- control$e

  # Compute the value of the objective (the negative Poisson
  # log-likelihood) at the initial iterate.
  f     <- cost(X,L,t(F),e,version = version)
  fbest <- f

  # Matrices in R are stored column-wise; to quickly access each row
  # of the matrix, the transpose of X is also stored.
  Xt <- t(X);
  
  # Compute the sum of the elements in each row and each column of the
  # counts matrix.
  xsrow <- rowSums(X)
  xscol <- colSums(X)
  
  # These are additional quantities used to implement the
  # extrapolation scheme. Here, "beta" and "betamax" are the
  # extrapolation parameters "beta" and "beta-bar" (the upper bound on
  # beta) used in Algorithm 3 of Ang & Gillis (2019).
  beta    <- 0
  betamax <- 0.99
  Fy      <- F
  Fn      <- F
  Fbest   <- F
  Ly      <- L
  Ln      <- L
  Lbest   <- L
  
  # Set up the data structure to record the algorithm's progress.
  progress <- as.matrix(data.frame(iter      = 1:numiter,
                                   objective = 0,
                                   mean.diff = 0,
                                   beta      = 0,
                                   timing    = 0))
  
  # Print a brief summary of the analysis, if requested.
  if (verbose) {
    cat(sprintf("Running %d EM + SQP updates ",numiter))
    cat("(fastTopics version 0.1-59)\n")
    if (control$extrapolate < Inf)
      cat(sprintf("Extrapolation begins at iteration %d\n",
                  control$extrapolate))
    else
      cat("Extrapolation is not active.\n")
    cat(sprintf("Data are %d x %d matrix with %0.1f%% nonzero proportion\n",
                n,m,100*mean(X > 0)))
    cat("iter objective (cost fn) mean.diff    beta\n")
  }

  # Iteratively apply the EM and SQP updates using the R or Rcpp
  # implementation.
  out <- altsqp_main_loop(X,Xt,F,Fn,Fy,Fbest,L,Ln,Ly,Lbest,f,fbest,xsrow,
                          xscol,beta,betamax,numiter,version,control,
                          progress,verbose)
  Fbest    <- out$Fbest
  Lbest    <- out$Lbest
  fbest    <- out$fbest
  progress <- out$progress
  rm(out)

  # Return a list containing (1) the estimate of the factors, (2) the
  # estimate of the loadings, (3) the value of the objective at these
  # estimates, and (4) a data frame recording the algorithm's progress
  # at each iteration.
  progress <- as.data.frame(progress)
  return(list(F = Fbest,L = Lbest,value = fbest,
              progress = as.data.frame(progress)))
}

# This helper function implements the main alt-SQP loop.
altsqp_main_loop <- function (X, Xt, F, Fn, Fy, Fbest, L, Ln, Ly, Lbest, f,
                              fbest, xsrow, xscol, beta, betamax, numiter,
                              version, control, progress, verbose) {
    
  # Get the optimization settings.
  nc               <- control$nc
  extrapolate      <- control$extrapolate
  beta.init        <- control$beta.init
  beta.increase    <- control$beta.increase
  beta.reduce      <- control$beta.reduce
  betamax.increase <- control$betamax.increase
  e                <- control$e
  
  # Iteratively apply the EM and SQP updates
  for (iter in 1:numiter) {

    # Store the value of the objective at the current iterate.
    f0 <- f

    # When the time is right, initiate the extrapolation scheme.
    if (beta == 0 & iter >= extrapolate) {
      beta  <- beta.init
      beta0 <- beta.init
    }

    timing <- system.time({

      # UPDATE FACTORS
      # --------------
      # Update the factors ("basis vectors").
      if (version == "Rcpp") {
        if (is.matrix(X))
          Fn <- t(altsqp_update_factors_rcpp(X,t(Fy),Ly,xscol,colSums(Ly),
                    e,control$numem,control$numsqp,control))
        else
          Fn <- t(altsqp_update_factors_sparse_rcpp(X,t(Fy),Ly,xscol,
                    colSums(Ly),e,control$numem,control$numsqp,control))
      } else {
        if (nc == 1)
          Fn <- altsqp.update.factors(X,Fy,Ly,xscol,control)
        else
          Fn <- altsqp.update.factors.multicore(X,Fy,Ly,xscol,control)
       }
      
      # Compute the extrapolated update for the factors. Note that
      # when beta = 0, Fy = Fn.
      Fy <- pmax(Fn + beta*(Fn - F),0)

      # UPDATE LOADINGS
      # ---------------
      # Update the loadings ("activations").
      if (version == "Rcpp") {
        if (is.matrix(Xt))
          Ln <- t(altsqp_update_loadings_rcpp(Xt,Fy,t(Ly),xsrow,colSums(Fy),
                    e,control$numem,control$numsqp,control))
        else
          Ln <- t(altsqp_update_loadings_sparse_rcpp(Xt,Fy,t(Ly),xsrow,
                    colSums(Fy),e,control$numem,control$numsqp,control))
      } else {
        if (nc == 1)
          Ln <- altsqp.update.loadings(Xt,Fy,Ly,xsrow,control)
        else
          Ln <- altsqp.update.loadings.multicore(Xt,Fy,Ly,xsrow,control)
      }
      
      # Compute the extrapolated update for the loadings. Note that
      # when beta = 0, Ly = Ln.
      Ly <- pmax(Ln + beta*(Ln - L),0)
    })

    # Compute the value of the objective (cost) function at the
    # extrapolated solution for the factors (F) and the
    # non-extrapolated solution for the loadings (L).
    f <- cost(X,Ln,t(Fy),e,version = version)
    if (beta == 0) {

      # No extrapolation is used, so use the basic coordinate-wise
      # updates for the factors and loadings.
      F <- Fn
      L <- Ln
    } else {

      # Update the extrapolation parameters following Algorithm 3 of
      # Ang & Gillis (2019).
      if (f > f0) {

        # The solution did not improve, so restart the extrapolation
        # scheme.
        Fy      <- F
        Ly      <- L
        betamax <- beta0
        beta    <- beta.reduce*beta
      } else {
        
        # The solution is improved; retain the basic co-ordinate ascent
        # update as well.
        F       <- Fn
        L       <- Ln
        beta    <- min(betamax,beta.increase * beta)
        beta0   <- beta
        betamax <- min(0.99,betamax.increase * betamax)
      }
    }        

    # If the solution is improved, update the current best solution
    # using the extrapolated estimates of the factors (F) and the
    # non-extrapolated estimates of the loadings (L).
    if (f < fbest) {
      d     <- mean((tcrossprod(Lbest,Fbest) - tcrossprod(Ln,Fy))^2)
      Fbest <- Fy
      Lbest <- Ln 
      fbest <- f
    } else
      d <- 0

    # Record the algorithm's progress.
    progress[iter,"objective"] <- fbest
    progress[iter,"mean.diff"] <- d
    progress[iter,"beta"]      <- beta
    progress[iter,"timing"]    <- timing["elapsed"]
    if (verbose)
      cat(sprintf("%4d %+0.12e %0.3e %0.1e\n",iter,fbest,d,beta))
  }

  return(list(Fbest = Fbest,Lbest = Lbest,fbest = fbest,progress = progress))
}

#' @rdname altsqp
#' 
#' @export
#' 
altsqp_control_default <- function()
  c(mixsqp_control_default(),
    list(nc               = 1,
         numem            = 1,
         numsqp           = 1,
         extrapolate      = 10,
         beta.init        = 0.5,
         beta.increase    = 1.1,
         beta.reduce      = 0.75,
         betamax.increase = 1.05))

# Update all the factors with the loadings remaining fixed.
altsqp.update.factors <- function (X, F, L, xscol, control) {
  m  <- ncol(X)
  ls <- colSums(L)
  for (j in 1:m) {
    i <- which(X[,j] > 0)
    if (control$numem > 0)
      F[j,] <- altsqp.update.em(L[i,],X[i,j],ls,xscol[j],F[j,],
                                control$e,control$numem)
    if (control$numsqp > 0)
      F[j,] <- altsqp.update.sqp(L[i,],X[i,j],ls,xscol[j],F[j,],
                                 control$numsqp,control)
  }
  return(F)
}

# Update all the loadings with the factors remaining fixed. Note that
# input argument X is the the *transpose* of the matrix inputted to
# altsqp (an m x n matrix).
altsqp.update.loadings <- function (X, F, L, xsrow, control) {
  n  <- ncol(X)
  fs <- colSums(F)
  for (i in 1:n) {
    j <- which(X[,i] > 0)
    if (control$numem > 0)
      L[i,] <- altsqp.update.em(F[j,],X[j,i],fs,xsrow[i],L[i,],
                                control$e,control$numem)
    if (control$numsqp > 0)
      L[i,] <- altsqp.update.sqp(F[j,],X[j,i],fs,xsrow[i],L[i,],
                                 control$numsqp,control)
  }
  return(L)  
}

# This is the multithreaded version of altsqp.update.factors,
# implemented using mclapply from the parallel package.
altsqp.update.factors.multicore <- function (X, F, L, xscol, control) {
  m    <- ncol(X)
  nc   <- control$nc
  cols <- splitIndices(m,nc)
  F <- mclapply(cols,
         function (j) altsqp.update.factors(X[,j],F[j,],L,xscol[j],control),
           mc.set.seed = FALSE,mc.allow.recursive = FALSE,mc.cores = nc)
  F <- do.call(rbind,F)
  F[unlist(cols),] <- F
  return(F)
}

# This is the multithreaded version of altsqp.update.loadings,
# implementing using mclapply from the parallel package. Note that
# input argument X is the the *transpose* of the matrix inputted to
# altsqp (an m x n matrix).
altsqp.update.loadings.multicore <- function (X, F, L, xsrow, control) {
  n    <- ncol(X)
  nc   <- control$nc
  rows <- splitIndices(n,nc)
  L <- mclapply(rows,
         function (i) altsqp.update.loadings(X[,i],F,L[i,],xsrow[i],control),
           mc.set.seed = FALSE,mc.allow.recursive = FALSE,mc.cores = nc)
  L <- do.call(rbind,L)
  L[unlist(rows),] <- L
  return(L)
}

# Run one EM update for the alternating SQP method.
altsqp.update.em <- function (B, w, bs, ws, x, e, numiter) {
  y   <- ws/bs
  out <- mixem(scale.cols(B,y),w/ws,x/y,numiter,e)
  return(out$x*y)
}

# Run one SQP update for the alternating SQP method.
altsqp.update.sqp <- function (B, w, bs, ws, x, numiter, control) {
  y   <- ws/bs
  out <- mixsqp(scale.cols(B,y),w/ws,x/y,numiter,control)
  return(out$x*y)
}
