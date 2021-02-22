#' @title Differential Count Analysis with a Multinomial Topic Model
#'
#' @description Implements methods for analysis of differential count
#' analysis using a topic model. These methods are motivated by gene
#' expression studies, but could have other uses, such as identifying
#' \dQuote{key words} in topics derived from text documents. To
#' improve accuracy of the differential expression analysis, an
#' empirical Bayes method is used to \dQuote{stabilize} the
#' estimates. A special case of \dQuote{hard} topic assignments is
#' also implemented---that is, the mixture proportions are all zeros
#' and ones---which involves greatly simplified (and faster)
#' calculations. Use \code{diff_count_clusters} for this special case.
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
#' \code{diff_count_analysis}), the Poisson model will closely
#' approximate a binomial model of the count data, so that the Poisson
#' model parameters \eqn{f_0, f_1} represent binomial probabilities.
#' In this case, \eqn{\beta} represents the LFC in \emph{relative}
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
#'   counts matrix to stabilize calculation of the LFC estimates and
#'   other statistics.
#'
#' @param fit.method Method used to estimate LFC and compute test
#'   statistics. The \code{"glm"} and \code{"optim"} computations are
#'   particularly slow, and it is recommended to use \code{fit.method =
#'   "em"} for most data sets. See \dQuote{Details} for more
#'   information.
#' 
#' @param shrink.method Method used to stabilize the LFC estimates.
#'   When \code{shrink.method = "ash"}, the "adaptive shrinkage" method
#'   implemented in the ashr package is used. When \code{shrink.method =
#'   "none"}, no stabilization is performed, and the \dQuote{raw} LFC
#'   estimates are returned.
#' 
#' @param numiter Maximum number of iterations performed in
#'   optimization of the Poisson model parameters. When
#'   \code{fit.method = "glm"}, this is passed as argument \code{maxit}
#'   to the \code{glm} function; when \code{fit.method = "optim"}, this
#'   is passed as argument \code{maxit} to the \code{optim} function.
#'
#' @param tol Controls the convergence tolerance for the optimization
#'   of the Poisson model parameters. When \code{fit.method = "glm"},
#'   this is passed as argument \code{epsilon} to function \code{glm};
#'   when \code{fit.method = "optim"}, the \code{factr} optimization
#'   setting is set to \code{tol * .Machine$double.eps}. When
#'   \code{fit.method = "em"}, the EM algorithm will be stopped early
#'   when the largest change between two successive updates is less than
#'   \code{tol}.
#' 
#' @param e A small, positive scalar included in some computations to
#'   avoid logarithms of zero and division by zero.
#'
#' @param show.warning Set \code{show.warning = FALSE} to suppress a
#'   message about calculations when mixture proportions are all 0 or 1.
#' 
#' @param verbose When \code{verbose = TRUE}, progress information is
#'   printed to the console.
#'
#' @param \dots When \code{shrink.method = "ash"}, these are
#'   additional arguments passed to \code{\link[ashr]{ash}}.
#' 
#' @return The return value is a list of m x k matrices, where m is
#'   the number of columns in the counts matrix, and k is the number of
#'   topics (for \code{diff_count_clusters}, m is the number of
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
#' @importFrom Matrix colMeans
#' @importFrom ashr ash
#' 
#' @export
#' 
diff_count_analysis <- function (fit, X, s = rowSums(X), pseudocount = 0.01,
                                 fit.method = c("em","optim","glm"),
                                 shrink.method = c("ash","none"),
                                 numiter = 100, tol = 1e-8, e = 1e-15,
                                 show.warning = TRUE, verbose = TRUE, ...) {

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
    stop("Input argument \"s\" should be a vector of positive numbers, with ",
         "length(s) = nrow(X)")
  normalize.by.totalcounts <- all(s == rowSums(X))
  
  # Check input argument "pseudocount".
  if (any(pseudocount <= 0))
    stop("Input argument \"pseudocount\" should be a positive number")

  # Process input arguments "fit.method" and "shrink.method".
  fit.method    <- match.arg(fit.method)
  shrink.method <- match.arg(shrink.method)
  
  # Get the number of topics (k) and the number of rows in the counts
  # matrix (m).
  m <- nrow(fit$F)
  k <- ncol(fit$F)

  # Compute the per-column averages. (This needs to be done before
  # we add the pseudocounts to the data.)
  colmeans <- colMeans(X)

  # Add the pseudocounts to the data.
  X <- rbind(X,matrix(pseudocount,k,ncol(X)))
  s <- c(s,rep(1,k))
  L <- rbind(fit$L,diag(k))

  # COMPUTE LOG-FOLD CHANGE STATISTICS
  # ----------------------------------
  # Fit the univariate ("single-count") Poisson models.
  if (verbose)
    cat(sprintf("Fitting %d x %d = %d univariate Poisson models.\n",m,k,m*k))

  # If all the mixture proportions are zeros or ones (or very close to
  # being zero or one), we can use the faster (analytical)
  # calculations implemented in fit_univar_poisson_models_hard;
  # otherwise, we call fit_univar_poisson_models.
  if (fit.method != "glm" & max(pmin(L,1-L)) <= 1e-14) {
    if (show.warning)
      warning("All mixture proportions are either zero or one; using ",
              "simpler single-topic calculations for model parameter ",
              "estimates")
    out <- fit_univar_poisson_models_hard(X,L,s,e)
  } else
    out <- fit_univar_poisson_models(X,L,s,
                                     ifelse(fit.method == "em","em-rcpp",
                                            fit.method),
                                     e,numiter,tol,verbose)
  if (normalize.by.totalcounts)
    if (any(out$F0 > 0.9) | any(out$F1 > 0.9))
      warning("One or more F0, F1 estimates are close to or greater than ",
              "1, so Poisson approximation to binomial may be poor; see ",
              "help(diff_count_analysis) for more information")
    
  # If glm was not used, compute the test statistics.
  if (fit.method != "glm") {
    if (verbose)
      cat("Computing log-fold change statistics.\n")
    out <- compute_univar_poisson_zscores_fast(X,L,out$F0,out$F1,s,e)
  }

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
  class(out) <- c("topic_model_diff_count","list")
  return(out)
}

#' @rdname diff_count_analysis
#'
#' @param cluster A factor, or vector that can be converted to a
#'   factor such as an integer or character vector, giving an assignment
#'   of samples (rows of \code{X}) to clusters. This could be, for
#'   example, the "cluster" output from \code{\link[stats]{kmeans}}.
#'
#' @param \dots Additional arguments passed to
#'   \code{diff_count_analysis}.
#' 
#' @export
#'
diff_count_clusters <- function (cluster, X, ...) {
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
  return(diff_count_analysis(fit,X,show.warning = FALSE,...))
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
