#' @title Differential Expression Analysis using a Topic Model
#'
#' @description Implements methods for analysis of differential count
#' analysis using a topic model. These methods are motivated by gene
#' expression studies, but could have other uses, such as identifying
#' \dQuote{key words} in topics derived from text documents. To
#' improve accuracy of the differential expression analysis, an
#' empirical Bayes method is used to \dQuote{stabilize} the
#' estimates. A special case of \dQuote{hard} topic assignments is
#' also implemented---that is, the topic proportions are all zeros
#' and ones---which involves greatly simplified (and faster)
#' calculations. Use \code{de_clusters} for this special case.
#'
#' @details The methods are based on the following univariate
#' (\dQuote{single-count}) Poisson model: \deqn{x_i ~ Poisson(s_i
#' \lambda_i),} in which the Poisson rates are defined as
#' \deqn{\lambda_i = (1 - q_i) f_0 + q_i f_1,} and where the two
#' unknowns to be estimated are \eqn{f_0, f_1 > 0.} This model is
#' applied separately to each count (column of the counts matrix), and
#' to each topic k. The \eqn{q_i}'s are the topic probabilities for a
#' given topic. An EM algorithm is used to compute maximum-likelihood
#' estimates (MLEs) of the two unknowns.
#'
#' The log-fold change (LFC) statistics are defined as the log-ratio
#' of the two Poisson rate parameters, \eqn{\beta = \log_2(f_1/f_0)}.
#' This statistic measures the increase (or decrease) in occurance in
#' one topic compared to all other topics. The use of the base-2
#' logarithm comes from the convention used in gene expression
#' studies.
#'
#' When each \eqn{s_i} (specified by input argument \code{s}) is equal
#' the total count for sample i (this is the default setting in
#' \code{de_analysis}), the Poisson model will closely approximate a
#' binomial model of the count data, so that the Poisson model
#' parameters \eqn{f_0, f_1} represent binomial probabilities.  In
#' this case, \eqn{\beta} represents the LFC in \emph{relative}
#' occurrence. (This Poisson approximation to the binomial is most
#' accurate when the total counts \code{rowSums(X)} are large and
#' \eqn{f_0, f_1} are small.)
#'
#' Other choices for \code{s} are possible, and implement different
#' normalization schemes for the counts, and different interpretations
#' of the LFC statistics. Setting \code{s} to all ones implements the
#' differential count analysis with no normalization, and \eqn{\beta}
#' represents the LFC in \emph{absolute} occurrence. This choice could
#' be appropriate in settings where the total count is well-controlled
#' across samples (rows of \code{X}).
#'
#' When \code{fit.method = "glm"}, the LFC estimates and test
#' statistics are implemented by \code{\link[stats]{glm}} with
#' \code{formula = x ~ b0 + b - 1} and \code{family = poisson(link =
#' "identity")}, where the unknowns are defined as \eqn{b_0 = f_0} and
#' \eqn{b = f_1 - f_0}. The outputted z-scores are p-values are for
#' the \code{b} coefficient; that is, they are test statistics for the
#' test that \eqn{f_0} and \eqn{f_1} are not equal.
#'
#' When \code{fit.method = "optim"} or \code{fit.method = "em"},
#' maximum-likelihood estimates are computed using
#' \code{\link[stats]{optim}} or using a fast EM algorithm, and test
#' statistics (standard errors, z-scores and p-values) are based on a
#' Laplace approximation to the likelihood at the MLE. These
#' calculations closely reproduce the test statistic calculations in
#' \code{\link[stats]{glm}}.
#'
#' We recommend setting \code{shrink.method = "ash"}, which uses the
#' \dQuote{adaptive shrinkage} method (Stephens, 2016) to improve
#' accuracy of the LFC estimates. The improvement in accuracy is
#' greatest for genes with low expression, or when the sample size is
#' small (or both). We follow the settings used in \code{lfcShrink}
#' from the DESeq2 package, with \code{type = "ashr"}.
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
#'   scaled in the Poisson model. See \dQuote{Details} for guidance on
#'   the choice of \code{s}.
#' 
#' @param pseudocount Observations with this value are added to the
#'   counts matrix to stabilize maximum-likelihood estimation.
#'
#' @param fit.method Method used to fit the Poisson models.  Note that
#'   \code{fit.method = "glm"} is the slowest method, and is mainly used
#'   for testing.
#'
#' @param lfc.stat Description of "lfc.stat" input argument goes here.
#' 
#' @param shrink.method Method used to stabilize the LFC estimates.
#'   When \code{shrink.method = "ash"}, the "adaptive shrinkage" method
#'   implemented in the ashr package is used. When \code{shrink.method =
#'   "none"}, no stabilization is performed, and the \dQuote{raw} LFC
#'   estimates are returned.
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
#' @param ns Describe input argument "ns" here.
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
#' 
#' @export
#' 
de_analysis <- function (fit, X, s = rowSums(X), pseudocount = 0.01,
                         fit.method = c("scd","em","mu","ccd","glm"),
                         lfc.stat = "le", shrink.method = c("ash","none"),
                         numiter = 20, minval = 1e-8, tol = 1e-8,
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
  
  # Ensure that none of the topic proportions are exactly zero or
  # exactly one.
  L <- fit$L
  L <- pmax(L,minval)
  L <- pmin(L,1 - minval)

  # SET UP DATA FOR POISSON MODELS
  # ------------------------------
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
  F <- fit_poisson_models(X,L,fit.method,eps,numiter,tol,nc)
  dimnames(F) <- dimnames(fit$F)

  # COMPUTE LOG-FOLD CHANGE STATISTICS
  # ----------------------------------
  # TO DO: Add some extra description here.
  if (verbose) {
    cat("Computing log-fold change statistics from ")
    cat(sprintf("%d Poisson models with k=%d.\n",m,k))
  }
  i <- c(5018,13171,13978,15685,sample(m,20))
  X <- X[,i]
  F <- F[i,]
  f0 <- f0[i]
  out <- compute_lfc_stats(X,F,L,f0,lfc.stat,ns,conf.level,rw,eps,nc,verbose)

  # Return the Poisson model MLEs (F), the log-fold change statistics
  # (est, low, high, z, lpval), and the relative rates under the "null"
  # model (f0).
  out$F <- F
  out$f0 <- f0
  class(out) <- c("topic_model_de_analysis","list")
  return(out)
  
  # STABILIZE ESTIMATES USING ADAPTIVE SHRINKAGE
  # --------------------------------------------
  # If requested, use "adaptive shrinkage" to stabilize the log-fold
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

  # PREPARE OUTPUTS
  # ---------------
  # Copy the row and column names used in "fit".
  rownames(out$F0)   <- rownames(fit$F)
  rownames(out$F1)   <- rownames(fit$F)
  rownames(out$beta) <- rownames(fit$F)
  rownames(out$se)   <- rownames(fit$F)
  rownames(out$Z)    <- rownames(fit$F)
  rownames(out$pval) <- rownames(fit$F)
  colnames(out$F0)   <- colnames(fit$F)
  colnames(out$F1)   <- colnames(fit$F)
  colnames(out$beta) <- colnames(fit$F)
  colnames(out$se)   <- colnames(fit$F)
  colnames(out$Z)    <- colnames(fit$F)
  colnames(out$pval) <- colnames(fit$F)

  # Return the Poisson model MLEs and the log-fold change statistics.
  out$colmeans <- colmeans
  class(out) <- c("topic_model_de_analysis","list")
  return(out)
}

#' @rdname de_analysis
#'
#' @param cluster A factor, or vector that can be converted to a
#'   factor such as an integer or character vector, giving an assignment
#'   of samples (rows of \code{X}) to clusters. This could be, for
#'   example, the "cluster" output from \code{\link[stats]{kmeans}}.
#'
#' @param \dots Additional arguments passed to \code{de_analysis}.
#' 
#' @export
#'
de_clusters <- function (cluster, X, ...) {
  if (!is.factor(cluster))
    cluster <- factor(cluster)

  # If necessary, remove all-zero columns from the counts matrix
  # before performing the calculations.
  if (any.allzero.cols(X)) {
    X <- remove.allzero.cols(X)
    warning(sprintf(paste("One or more columns of X are all zero; after",
                          "removing all-zero columns, %d columns will be",
                          "used for all statistical calculations"),ncol(X)))
  }

  fit <- fit_multinom_model(cluster,X)
  return(de_analysis(fit,X,...))
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

  # Output the posterior estimates (b), the posterior standard errors
  # (se), the z-scores (Z) and -log10 p-values (pval).
  return(out)
}
