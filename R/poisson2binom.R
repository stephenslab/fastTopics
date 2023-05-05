#' @title Recover Binomial Topic Model Fit from Poisson NMF fit
#'
#' @description Add a brief description
#'
#' @details Describe the binomial topic model in more detail here.
#'  
#' @param X The n x m \dQuote{binary} matrix; all entries of X should
#'   be between 0 and 1 (including 0 and 1). It can be a sparse matrix
#'   (class \code{"dgCMatrix"}) or dense matrix (class \code{"matrix"}).
#'
#' @param fit Describe input argument "fit" here.
#'
#' @param numem Describe input argument "numem" here.
#'
#' @param umin Describe input argument "umin" here.
#' 
#' @param verbose Describe input argument "verbose" here.
#' 
#' @return Describe the output here.
#' 
#' @examples
#' # See the vignette for an example.
#' 
#' @export
#'
poisson2binom <- function (X, fit, numem = 0, umin = 1e-4, verbose = TRUE) {
  if (!requireNamespace("NNLM",quietly = TRUE))
    stop("poisson2binom requires the NNLM package")

  # Check input argument "fit".
  if (inherits(fit,"binom_topic_model_fit"))
    return(fit)
  if (!inherits(fit,"poisson_nmf_fit"))
    stop("Input argument \"fit\" should be an object of class ",
         "\"poisson_nmf_fit\"")
  verify.fit(fit)
  if (ncol(fit$F) < 2 | ncol(fit$L) < 2)
    stop("Input matrices \"fit$F\" and \"fit$L\" should have 2 or more",
         "columns")

  # Check and process input argument "X".
  if (!((is.numeric(X) & is.matrix(X)) | is.sparse.matrix(X)))
    stop("Input argument \"X\" should be a numeric matrix (a \"matrix\" or ",
         "a \"dgCMatrix\")")
  verify.fit.and.count.matrix(X,fit)
  if (any(X < 0) | any(X > 1))
    warning("Input argument \"X\" should be a \"binary\" matrix ",
            "(that is, all entries should range from 0 and 1)")
  if (is.matrix(X) & is.integer(X))
    storage.mode(X) <- "double"

  # Choose U = diag(u) such that L*U is closer to being a matrix of
  # topic proportions.
  ones <- matrix(1,n,1)
  L    <- fit$L
  F    <- fit$F
  if (verbose)
    cat("Rescaling L and F using non-negative linear regression (nnlm).\n")
  u   <- drop(coef(NNLM::nnlm(L,ones)))
  u   <- pmax(u,umin)
  # TO DO: Make sure u's are always positive.
  L   <- scale.cols(L,u)
  L   <- normalize.rows(L)
  F   <- scale.cols(F,1/u)
  fit <- list(F = F,L = L,progress = NA)

  # Refine the binomial topic model fit by performing several EM updates.
  if (numem > 0) {
    cat("Performing",numem,"EM updates to refine the fit.\n")
    progress <- as.data.frame(cbind(1:numem,0,0))
    names(progress) <- c("iter","delta.f","delta.l")
    if (verbose)
    cat("iter  |F - F'|  |L - L'|\n")
    for (i in 1:numem) {
      fit0 <- fit
      fit  <- fit_binom_topic_model_em(X,fit,numem)
      progress[i,"delta.f"] <- max(abs(fit0$F - fit$F))
      progress[i,"delta.l"] <- max(abs(fit0$L - fit$L))
      if (verbose)
        cat(sprintf("%4d %0.3e %0.3e\n",i,progress[i,"delta.f"],
                    progress[i,"delta.f"]))
    }
    fit$progress <- progress
  }
  
  # Return the Binomial topic model fit.
  class(fit) <- c("binom_topic_model_fit","list")
  return(fit)
}

# Perform a single EM udpate for fiitting the binomial topic model to
# binary data matrix X. This code is adapted from the meth_tpxEM
# function in the methClust package by Kushal Dey.
fit_binom_topic_model_em <- function (X, fit, numiter) {
  if (!is.matrix(X))
    X <- as.matrix(X)
    
  # Make sure no parameters are exactly zero or exactly one.
  e <- 1e-8
  L <- fit$L
  F <- fit$F
  F <- clamp(F,e,1 - e)
  L <- clamp(L,e,1 - e)
  L <- normalize.rows(L)

  # Perform the E step.
  A  <- X/tcrossprod(L,F)
  M  <- (A %*% F) * L
  Mt <- crossprod(A,L) * F
  A  <- (1 - X)/tcrossprod(L,1 - F)
  U  <- (A %*% (1 - F)) * L
  Ut <- crossprod(A,L) * (1 - F)

  # Perform the M step.
  L <- normalize.rows(M + U)
  F <- Mt/(Mt + Ut)
  return(list(F = F,L = L))
}
