#' @title Differential Expression Analysis using a Topic Model
#'
#' @description Implements methods for differential expression
#'   analysis using a topic model. These methods are motivated by gene
#'   expression studies, but could have other uses, such as identifying
#'   \dQuote{key words} for topics.
#'
#' @details The methods are based on the Poisson model \deqn{x_i ~
#' Poisson(s_i \lambda_i),} in which the Poisson rates are
#' \deqn{\lambda_i = \sum_{j=1}^k s_i l_{ij} f_j,} in which the
#' \eqn{l_{ik}} are the topic proportions and the \eqn{f_j} are the
#' unknowns to be estimated. This model is applied separately to each
#' column of \code{X}. When \eqn{s_i} (specified by input argument
#' \code{s}) is equal the total count in row i (this is the default),
#' the Poisson model will closely approximate a binomial model of the
#' count data, and the unknowns \eqn{f_j} will approximate binomial
#' probabilities. (The Poisson approximation to the binomial is most
#' accurate when the total counts \code{rowSums(X)} are large and the
#' unknowns \eqn{f_j} are small.) Other choices for \code{s} are
#' possible, and implement different normalization schemes.
#'
#' TO DO: Describe LFC statistics here.
#' 
#' We recommend setting \code{shrink.method = "ash"}, which uses the
#' \dQuote{adaptive shrinkage} method (Stephens, 2016) to improve
#' accuracy of the LFC estimates. We follow the settings used in
#' \code{lfcShrink} from the DESeq2 package, with \code{type =
#' "ashr"}.
#' 
#' @param fit An object of class \dQuote{poisson_nmf_fit} or
#'   \dQuote{multinom_topic_model_fit}. If a Poisson NMF fit is provided
#'   as input, the corresponding multinomial topic model fit is
#'   automatically recovered using \code{\link{poisson2multinom}}.
#'
#' @param X The n x m counts matrix. It can be a sparse matrix (class
#'   \code{"dgCMatrix"}) or dense matrix (class \code{"matrix"}).
#'
#' @param s A numeric vector of length n determining how the rates are
#'   scaled in the Poisson models. See \dQuote{Details} for guidance on
#'   the choice of \code{s}.
#' 
#' @param pseudocount Observations with this value are added to the
#'   counts matrix to stabilize maximum-likelihood estimation.
#'
#' @param fit.method Method used to fit the Poisson models. Note that
#'   \code{fit.method = "glm"} is the slowest method, and is mainly used
#'   for testing.
#'
#' @param shrink.method Method used to stabilize the LFC estimates.
#'   When \code{shrink.method = "ash"}, the "adaptive shrinkage" method
#'   implemented in the ashr package is used. When \code{shrink.method =
#'   "none"}, no stabilization is performed, and the \dQuote{raw} LFC
#'   estimates are returned.
#'
#' @param lfc.stat The log-fold change statistics returned:
#'   \code{lfc.stat = "vsnull"}, the log-fold change relative to the
#'   null; \code{lfc.stat = le"}, the \dQuote{least extreme} LFC; or a
#'   topic name or number, in which case the LFC is defined relative to
#'   the selected topic. See \dQuote{Details} for more detailed
#'   explanations of these choices.
#' 
#' @param numiter Maximum number of iterations performed in fitting
#'   the Poisson models. When \code{fit.method = "glm"}, this is passed
#'   as argument \code{maxit} to the \code{glm} function.
#'
#' @param minval A small, positive number. All topic proportions less
#'   than this value and greater than \code{1-minval} are set to this
#'   value.
#' 
#' @param tol Controls the convergence tolerance for fitting the
#'   Poisson models. When \code{fit.method = "glm"}, this is passed as
#'   argument \code{epsilon} to function \code{glm}.
#'
#' @param conf.level Describe input argument "conf.level" here.
#' 
#' @param ns Number of Monte Carlo samples simulated by random-walk
#'   MCMC for estimating posterior LFC quantities.
#'
#' @param rw Describe input argument "rw" here.
#' 
#' @param eps A small, non-negative number added to the terms inside
#'   the logarithms to avoid computing logarithms of zero.
#'
#' @param nc Number of threads used in the multithreaded computations.
#' 
#' @param verbose When \code{verbose = TRUE}, progress information is
#'   printed to the console.
#'
#' @param \dots When \code{shrink.method = "ash"}, these are
#'   additional arguments passed to \code{\link[ashr]{ash}}.
#'
#' @return The return value is a list of m x k matrices, where m is
#'   the number of columns in the counts matrix, and k is the number of
#'   topics (for \code{de_clusters}, m is the number of
#'   clusters), and an additional vector:
#'
#' \item{colmeans}{A vector of length m containing the count averages
#'   (\code{colMeans(X)}).}
#'
#' \item{F0}{Estimates of the Poisson model parameters \eqn{f_0}.}
#'
#' \item{F1}{Estimates of the Poisson model parameters \eqn{f_1}.}
#'
#' \item{beta}{LFC estimates \code{beta = log2(F1/F0)}.}
#'
#' \item{se}{Standard errors of the Poisson glm parameters \eqn{b =
#'   f_1 - f_0}. }
#'
#' \item{Z}{z-scores for the Poisson glm parameters \eqn{b = f_1 -
#'   f_0}.}
#'
#' \item{pval}{-log10 two-tailed p-values computed from the
#'   z-scores. In some extreme cases the calculations may produce zero
#'   or negative standard errors, in which case the standard errors are
#'   set to \code{NA} and the z-scores and -log10 p-values are set to
#'   zero.}
#'
#' @references
#' Stephens, M. (2016). False discovery rates: a new deal.
#' \emph{Biostatistics} \bold{18}(2), kxw041.
#' \url{https://doi.org/10.1093/biostatistics/kxw041}
#'
#' Zhu, A., Ibrahim, J. G. and Love, M. I. (2019). Heavy-tailed prior
#' distributions for sequence count data: removing the noise and
#' preserving large differences. \emph{Bioinformatics} \bold{35}(12),
#' 2084–2092.
#' 
#' @importFrom Matrix rowSums
#' @importFrom Matrix colSums
#' @importFrom ashr ash
#' @importFrom stats rnorm
#' @importFrom stats runif
#'
#' @export
#' 
de_analysis <- function (fit, X, s = rowSums(X), pseudocount = 0.01,
                         fit.method = c("scd","em","mu","ccd","glm"),
                         shrink.method = c("ash","none"), lfc.stat = "de",
                         numiter = 20, minval = 1e-10, tol = 1e-8,
                         conf.level = 0.9, ns = 1000, rw = 0.3, eps = 1e-15,
                         nc = 1, verbose = TRUE, ...) {

  # CHECK AND PROCESS INPUTS
  # ------------------------
  # Check and process input argument "fit".
  if (!(inherits(fit,"poisson_nmf_fit") |
        inherits(fit,"multinom_topic_model_fit")))
    stop("Input \"fit\" should be an object of class \"poisson_nmf_fit\" or ",
         "\"multinom_topic_model_fit\"")
  if (inherits(fit,"poisson_nmf_fit"))
    fit <- poisson2multinom(fit)
  
  # Check and process input argument "X".
  if (!((is.numeric(X) & is.matrix(X)) | is.sparse.matrix(X)))
    stop("Input argument \"X\" should be a numeric matrix (a \"matrix\" or ",
         "a \"dgCMatrix\")")
  verify.fit.and.count.matrix(X,fit)
  if (is.matrix(X) & is.integer(X))
    storage.mode(X) <- "double"

  # Check input argument "s".
  verify.positive.vector(s)
  if (length(s) != nrow(X))
    stop("Input argument \"s\" should be a vector of positive numbers, ",
         "in which length(s) = nrow(X)")
  
  # Check input argument "pseudocount".
  if (any(pseudocount <= 0))
    stop("Input argument \"pseudocount\" should be a positive number")

  # Process input arguments "fit.method" and "shrink.method".
  fit.method <- match.arg(fit.method)
  shrink.method <- match.arg(shrink.method)

  # Process input argument "lfc.stat".
  # TO DO.
  
  # Get the number of rows (n) and columns (m) in the counts matrix, and
  # the number of groups or topics (k).
  n <- nrow(X)
  m <- ncol(X)
  k <- ncol(fit$F)

  # FIT NULL MODELS
  # ---------------
  # Compute the MLE f0 in the "null" model x ~ Poisson(u), with u =
  # s*f0. (This calculation is performed for each column of the counts
  # matrix. It must be done before adding pseudocounts to the data.)
  # Equivalently, this is the maximum-likelihood estimate of the
  # binomial probability in the "null" Binomial model x ~ Binom(s,p0).
  f0 <- colSums(X)/sum(s)
  names(f0) <- rownames(fit$F)
  
  # SET UP DATA FOR FITTING POISSON MODELS
  # --------------------------------------
  # Ensure that none of the topic proportions are exactly zero or
  # exactly one.
  L <- fit$L
  L <- pmax(L,minval)
  L <- pmin(L,1 - minval)

  # Add "pseudocounts" to the data, and get the Poisson NMF loadings
  # matrix. From this point on, we will fit Poisson glm models x ~
  # Poisson (u), u = sum(L*f), where x is a column of X.
  out <- add_pseudocounts(X,s*L,pseudocount)
  X <- out$X
  L <- out$L

  # FIT POISSON MODELS
  # ------------------
  # For each column j of the counts matrix, compute MLEs of the
  # parameters in the Poisson glm, x ~ Poisson(u), in which the
  # Poisson rates are u = sum(L*f), and f = F[j,].
  if (verbose)
    cat(sprintf("Fitting %d Poisson models with k=%d using method=\"%s\".\n",
                m,k,fit.method))
  nc <- initialize.multithreading(nc)
  F <- fit_poisson_models(X,L,fit.method,eps,numiter,tol,nc)
  F <- pmax(F,minval)
  dimnames(F) <- dimnames(fit$F)

  # COMPUTE LOG-FOLD CHANGE STATISTICS
  # ----------------------------------
  # Perform MCMC to simulate the posterior distribution of the LFC
  # statistics, then compute key posterior quantities from the
  # simulated Monte Carlo samples.
  if (verbose) {
    cat("Computing log-fold change statistics from ")
    cat(sprintf("%d Poisson models with k=%d.\n",m,k))
  }
  D <- matrix(rnorm(ns*k),ns,k)
  U <- matrix(runif(ns*k),ns,k)
  if (nc == 1)
    out <- compute_lfc_stats(X,F,L,f0,D,U,lfc.stat,conf.level,rw,eps,verbose)
  else {
    message(sprintf("Using %d SOCK threads.",nc))
    if (verbose)
      message("Progress bar cannot be shown when multithreading is used. ",
              "For large data sets, the total runtime may be estimated by ",
              "performing an initial test run with a small number of Monte ",
              "Carlo samples, e.g., ns = 10; the total runtime should scale ",
              "roughly linearly in ns, the number of Monte Carlo samples.")
    out <- compute_lfc_stats_multicore(X,F,L,f0,D,U,lfc.stat,conf.level,
                                       rw,eps,nc)
  }

  # STABILIZE ESTIMATES USING ADAPTIVE SHRINKAGE
  # --------------------------------------------
  # If requested, use adaptive shrinkage to stabilize the log-fold
  # change estimates.
  if (shrink.method == "ash") {
    if (verbose)
      cat("Stabilizing log-fold change estimates using adaptive shrinkage.\n")
    res <- shrink_lfc(out$F1 - out$F0,out$se,...)
    out$F1   <- out$F0 + res$b
    out$beta <- log2(out$F1/out$F0)
    out$se   <- res$se
    out$Z    <- res$Z
    out$pval <- res$pval
  }
  
  # Return the Poisson model MLEs (F), the log-fold change statistics
  # (est, low, high, z, lpval), and the relative rates under the "null"
  # model (f0).
  out$F <- F
  out$f0 <- f0
  class(out) <- c("topic_model_de_analysis","list")
  return(out)
}

# Perform adaptive shrinkage on the unknowns b = f1 - f0.
shrink_lfc <- function (b, se, ...) {

  # Get the number of effect estimates (m) and the number of topics (k).
  m <- nrow(b)
  k <- ncol(b)

  # Initialize the outputs.
  out <- list(b    = b,
              se   = se,
              Z    = matrix(0,m,k),
              pval = matrix(0,m,k))

  # Repeat for each topic.
  for (j in 1:k) {
    i <- which(!is.na(out$se[,j]))
    if (length(i) > 0) {

      # Run adaptive shrinkage, then extract the posterior estimates (b)
      # and standard errors (se).
      ans <- ash(out$b[i,j],out$se[i,j],mixcompdist = "normal",
                 method = "shrink",...)
      b   <- ans$result$PosteriorMean
      se  <- ans$result$PosteriorSD
      z   <- b/se

      # Store the stabilized estimates.
      out$b[i,j]    <- b
      out$se[i,j]   <- se
      out$Z[i,j]    <- z
      out$pval[i,j] <- -lpfromz(z)
    }
  }

  # TO DO: Add row and column names to the outputs, if necessary.
  
  # Output the posterior estimates (b), the posterior standard errors
  # (se), the z-scores (Z) and -log10 p-values (pval).
  return(out)
}
