#' @title Alternating SQP Method for Optimizing Topic Models and Non-negative Matrix Factorizations
#' 
#' @description Compute maximum-likelihood estimates for the Poisson
#'   topic model with factors F and loadings L; equivalently, find a
#'   non-negative matrix factorization \eqn{L * F^T} of matrix X that
#'   optimizes the beta (or Bregman) divergence objective.
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
#' \item{\code{nc}}{Describe nc here.}
#'
#' \item{\code{extrapolate}}{Describe extrapolate here.}
#'
#' \item{\code{beta0}}{Describe beta0 here.}
#'
#' \item{\code{betamaxinc}}{Describe betamaxinc here.}
#'
#' \item{\code{betainc}}{Describe betainc here.}
#'
#' \item{\code{betared}}{Describe betared here.}
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
#'   small inaccuracy in the computation.
#'
#' @return A list object with the following elements:
#'
#' \item{F}{Describe F here.}
#'
#' \item{L}{Describe L here.}
#'
#' \item{value}{Describe "value" here.}
#'
#' \item{progress}{Describe "progress" here.}
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
#' X <- matrix(rpois(n*m,L %*% t(F)),n,m)
#' X <- as(X,"sparseMatrix")
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
#' fit2 <- altsqp(X,fit0,numiter = 60,verbose = FALSE)
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
#' @importFrom utils modifyList
#' @importFrom Matrix rowSums
#' @importFrom Matrix colSums
#' @importFrom Matrix mean
#' @importFrom parallel splitIndices
#' @importFrom parallel mclapply
#' 
#' @export
#' 
altsqp <- function (X, fit, numiter = 100, control = list(), verbose = TRUE) {

  # Verify and process input matrix X. Each row and each column of the
  # matrix should have at least two positive entries.
  verify.matrix(X)
  if (!(all(rowSums(X > 0) >= 2) &
        all(colSums(X > 0) >= 2)))
    stop(paste("Each row and column of \"X\" should have at least two",
               "positive entries"))
  if (!inherits(X,"sparseMatrix") & mean(X > 0) < 0.1)
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

  # Verify and process input matrix L.
  verify.matrix(fit$L)
  if (!is.matrix(L))
    L <- as.matrix(L)
  if (is.integer(L))
    storage.mode(L) <- "double"

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
  
  # Get the optimization settings.
  control <- modifyList(altsqp_control_default(),control,keep.null = TRUE)
  nc               <- control$nc
  extrapolate      <- control$extrapolate
  beta.init        <- control$beta.init
  betamax.increase <- control$betamax.increase
  beta.increase    <- control$beta.increase
  beta.reduce      <- control$beta.reduce
  e                <- control$e

  # Compute the value of the objective (the negative Poisson
  # log-likelihood) at the initial iterate.
  f     <- cost(X,tcrossprod(L,F),e)
  fbest <- f
  
  # These are additional quantities used to implement the
  # extrapolation scheme. Here, "beta" and "betamax" are the
  # extrapolation parameters "beta" and "beta-bar" (the upper bound on
  # beta) used in Algorithm 3 of Ang & Gillis (2019).
  beta     <- 0
  betamax  <- 0.99
  Fy       <- F
  Fbest    <- F
  Ly       <- L 
  Lbest    <- L
  
  # Set up the data structure to record the algorithm's progress.
  progress <- data.frame(iter      = 1:numiter,
                         objective = 0,
                         max.diff  = 0,
                         beta      = 0,
                         timing    = 0)

  # Print a brief summary of the analysis, if requested.
  if (verbose) {
    cat(sprintf("Running %d alternating SQP updates ",numiter))
    cat("(fastTopics version 0.1-20)\n")
    if (extrapolate < Inf)
      cat(sprintf("Extrapolation begins at iteration %d\n",extrapolate))
    else
      cat("Extrapolation is not active.\n")
    cat(sprintf("Data are %d x %d matrix with %0.1f%% nonzero proportion\n",
                n,m,100*mean(X > 0)))
    cat("iter         objective max.diff    beta\n")
  }
  
  # Iteratively apply the EM And SQP updates.
  for (iter in 1:numiter) {

    # Store the value of the objective at the current iterate, and the
    # current setting of the extrapolation parameter.
    f0    <- f
    beta0 <- beta

    # When the time is right, initiate the extrapolation scheme.
    if (beta == 0 & iter >= extrapolate)
      beta <- beta.init

    timing <- system.time({

      # UPDATE FACTORS
      # --------------
      # Update the factors ("basis vectors").
      if (nc == 1)
        Fn <- altsqp.update.factors(X,Fy,Ly,control)
      else
        Fn <- altsqp.update.factors.multicore(X,Fy,Ly,control)

      # Compute the extrapolated update for the factors.
      Fy <- pmax(Fn + beta*(Fn - F),0)

      # UPDATE LOADINGS
      # ---------------
      # Update the loadings ("activations").
      if (nc == 1)
        Ln <- altsqp.update.loadings(X,Fy,Ly,control)
      else
        Ln <- altsqp.update.loadings.multicore(X,Fy,Ly,control)

      # Compute the extrapolated update for the loadings.
      Ly <- pmax(Ln + beta*(Ln - L),0)
    })

    # Compute the value of the objective (cost) function at the
    # extrapolated solution for the factors (F) and the
    # non-extrapolated solution for the loadings (L).
    f <- cost(X,tcrossprod(Ln,Fy),e)

    if (beta == 0) {
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
        betamax <- min(0.99,betamax.increase * betamax)
      }
    }        

    # If the solution is improved, update the current best solution.
    if (f < fbest) {
      d     <- max(abs(tcrossprod(Lbest,Fbest) - tcrossprod(Ln,Fy)))
      Fbest <- Fy
      Lbest <- Ln 
      fbest <- f
    } else
      d <- 0

    # Record the algorithm's progress.
    progress[iter,"objective"] <- fbest
    progress[iter,"max.diff"]  <- d
    progress[iter,"beta"]      <- beta
    progress[iter,"timing"]    <- timing["elapsed"]
    if (verbose)
      cat(sprintf("%4d %+0.10e %0.2e %0.1e\n",iter,fbest,d,beta))
  }

  # Return a list containing (1) the estimate of the factors, (2) the
  # estimate of the loadings, (3) the value of the objective at these
  # estimates, and (4) a data frame recording the algorithm's progress
  # at each iteration.
  return(list(F = Fbest,L = Lbest,value = fbest,progress = progress))
}

#' @rdname altsqp
#' 
#' @export
#' 
altsqp_control_default <- function()
  c(mixsqp_control_default(),
    list(nc = 1,extrapolate = 10,beta.init = 0.5,betamax.increase = 1.05,
         beta.increase = 1.1,beta.reduce = 0.75))

# Update all the factors with the loadings remaining fixed.
altsqp.update.factors <- function (X, F, L, control) {
  m <- ncol(X)
  for (j in 1:m) {
    F[j,] <- altsqp.update.em(L,X[,j],F[j,],control$e)
    F[j,] <- altsqp.update.sqp(L,X[,j],F[j,],control)
  }
  return(F)
}

# Update all the loadings with the factors remaining fixed.
altsqp.update.loadings <- function (X, F, L, control) {
  n <- nrow(X)
  for (i in 1:n) {
    L[i,] <- altsqp.update.em(F,X[i,],L[i,],control$e)
    L[i,] <- altsqp.update.sqp(F,X[i,],L[i,],control)
  }
  return(L)  
}

# This is the multithreaded version of altsqp.update.factors,
# implemented using mclapply from the parallel package.
altsqp.update.factors.multicore <- function (X, F, L, control) {
  cols <- splitIndices(m,nc)
  F <- mclapply(cols,
         function (j) altsqp.update.factors(X[,j],F[j,],L,control),
           mc.set.seed = FALSE,mc.allow.recursive = FALSE,mc.cores = nc)
  F <- do.call(rbind,F)
  F[unlist(cols),] <- F
  return(F)
}

# This is the multithreaded version of altsqp.update.loadings,
# implementing using mclapply from the parallel package.
altsqp.update.loadings.multicore <- function (X, F, L, control) {
  rows <- splitIndices(n,nc)
  L <- mclapply(rows,
         function (i) altsqp.update.loadings(X[i,],F,L[i,],control),
           mc.set.seed = FALSE,mc.allow.recursive = FALSE,mc.cores = nc)
  L <- do.call(rbind,L)
  L[unlist(rows),] <- L
  return(L)
}

# Run one EM update for the alternating SQP method.
altsqp.update.em <- function (B, w, y, e) {

  # Remove any counts that are exactly zero.
  ws <- sum(w)
  bs <- colSums(B)
  i  <- which(w > 0)
  w  <- w[i]
  B  <- B[i,]

  # -----
  # TO DO: Verify that all inputs to mixem are valid.
  # -----
  
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

  # -----
  # TO DO: Verify that all inputs to mixsqp are valid.
  # -----
  
  # Run an SQP update for the modified problem.
  out <- mixsqp(scale.cols(B,ws/bs),w/ws,y*bs/ws,1,control)

  # Recover the solution to the unmodified problem.
  return(out$x*ws/bs)
}
