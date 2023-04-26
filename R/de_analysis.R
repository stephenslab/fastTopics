#' @rdname de_analysis
#' 
#' @title Differential Expression Analysis using a Topic Model
#'
#' @description Implements methods for differential expression
#'   analysis using a topic model. These methods are motivated by gene
#'   expression studies, but could have other uses, such as identifying
#'   \dQuote{key words} for topics.
#'
#' @details The methods are based on the Poisson model
#' \deqn{x_i ~ Poisson(\lambda_i),} in which the Poisson rates are
#' \deqn{\lambda_i = \sum_{j=1}^k s_i l_{ij} f_j,} the \eqn{l_{ik}}
#' are the topic proportions and the \eqn{f_j} are the unknowns to be
#' estimated. This model is applied separately to each column of
#' \code{X}. When \eqn{s_i} (specified by input argument \code{s}) is
#' equal the total count in row i (this is the default), the Poisson
#' model will closely approximate a binomial model of the count data,
#' and the unknowns \eqn{f_j} will approximate binomial
#' probabilities. (The Poisson approximation to the binomial is most
#' accurate when the total counts \code{rowSums(X)} are large and the
#' unknowns \eqn{f_j} are small.) Other choices for \code{s} are
#' possible, and implement different normalization schemes.
#'
#' To allow for some flexibility, \code{de_analysis} allows for the
#' log-fold change to be measured in several ways.
#'
#' One option is to compare against the probability under the null
#' model: \eqn{LFC(j) = log2(f_j/f_0)}, where \eqn{f_0} is the single
#' parameter in the Poisson model \eqn{x_i ~ Poisson(\lambda_i)} with
#' rates \eqn{\lambda_i = s_i f_0}. This LFC definition is chosen with
#' \code{lfc.stat = "vsnull"}.
#'
#' Another option is to compare against a chosen topic, k: \eqn{LFC(j)
#' = log2(f_j/f_k)}. By definition, \eqn{LFC(k)} is zero, and
#' statistics such as z-scores and p-values for topic k are set to
#' \code{NA}. This LFC definition is selected by setting
#' \code{lfc.stat = k}.
#'
#' A final option (which is the default) computes the \dQuote{least
#' extreme} LFC, defined as \eqn{LFC(j) = log2(f_j/f_k)} such that
#' \eqn{k} is the topic other than \eqn{j} that gives the ratio
#' \eqn{f_j/f_k} closest to 1. This option is chosen with
#' \code{lfc.stat = "le"}.
#' 
#' We recommend setting \code{shrink.method = "ash"}, which uses the
#' \dQuote{adaptive shrinkage} method (Stephens, 2016) to improve
#' accuracy of the posterior mean estimates and z-scores. We follow
#' the settings used in \code{lfcShrink} from the \sQuote{DESeq2}
#' package, with \code{type = "ashr"}.
#'
#' Note that all LFC statistics are defined using the base-2 logarithm
#' following the conventioned used in differential expression
#' analysis.
#' 
#' The \code{control} argument is a list in which any of the
#' following named components will override the default optimization
#' algorithm settings (as they are defined by
#' \code{de_analysis_control_default}):
#' 
#' \describe{
#'
#' \item{\code{numiter}}{Maximum number of iterations performed in
#'   fitting the Poisson models. When \code{fit.method = "glm"}, this is
#'   passed as argument \code{maxit} to the \code{glm} function.}
#'
#' \item{\code{minval}}{A small, positive number. All topic
#'   proportions less than this value and greater than \code{1 - minval}
#'   are set to this value.}
#' 
#' \item{\code{tol}}{Controls the convergence tolerance for fitting
#'   the Poisson models. When \code{fit.method = "glm"}, this is passed
#'   as argument \code{epsilon} to function \code{glm}.}
#'
#' \item{\code{conf.level}}{The size of the highest posterior density
#'   (HPD) intervals. Should be a number greater than 0 and less than 1.}
#' 
#' \item{\code{ns}}{Number of Monte Carlo samples simulated by
#'   random-walk MCMC for estimating posterior LFC quantities.}
#'
#' \item{\code{rw}}{The standard deviation of the normal density used
#'   to propose new states in the random-walk MCMC.}
#' 
#' \item{\code{eps}}{A small, non-negative number added to the terms
#'   inside the logarithms to avoid computing logarithms of zero.}
#'
#' \item{\code{nc}}{Number of threads used in the multithreaded
#'   computations. Note that the multithreading relies on forking hence
#'   is not available on Windows; will return an error on Windows
#'   unless \code{nc = 1}. See \code{\link[parallel]{mclapply}} for
#'   details.}
#'
#' \item{\code{nsplit}}{The number of data splits used in the
#'   multithreaded computations (only relevant when \code{nc > 1}). More
#'   splits increase the granularity of the progress bar, but can also
#'   slow down the mutithreaded computations by introducing more
#'   overhead in the call to \code{\link[pbapply]{pblapply}}.}}
#' 
#' @param fit An object of class \dQuote{poisson_nmf_fit} or
#'   \dQuote{multinom_topic_model_fit}, or an n x k matrix of topic
#'   proportions, where k is the number of topics. (The elements in each
#'   row of this matrix should sum to 1.) If a Poisson NMF fit is
#'   provided as input, the corresponding multinomial topic model fit is
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
#' @param shrink.method Method used to stabilize the posterior mean
#'   LFC estimates.  When \code{shrink.method = "ash"}, the "adaptive
#'   shrinkage" method implemented in the \sQuote{ashr} package is used to
#'   compute posterior. When \code{shrink.method = "none"}, no
#'   stabilization is performed, and the \dQuote{raw} posterior mean LFC
#'   estimates are returned.
#'
#' @param lfc.stat The log-fold change statistics returned:
#'   \code{lfc.stat = "vsnull"}, the log-fold change relative to the
#'   null; \code{lfc.stat = "le"}, the \dQuote{least extreme} LFC; or a
#'   topic name or number, in which case the LFC is defined relative to
#'   the selected topic. See \dQuote{Details} for more detailed
#'   explanations of these choices.
#' 
#' @param control A list of parameters controlling behaviour of
#'   the optimization and Monte Carlo algorithms. See \sQuote{Details}.
#' 
#' @param verbose When \code{verbose = TRUE}, progress information is
#'   printed to the console.
#'
#' @param \dots When \code{shrink.method = "ash"}, these are
#'   additional arguments passed to \code{\link[ashr]{ash}}.
#'
#' @return A list with the following elements:
#'
#' \item{est}{The log-fold change maximum-likelihood estimates.}
#'
#' \item{postmean}{Posterior mean LFC estimates.}
#' 
#' \item{lower}{Lower limits of estimated HPD intervals. Note that
#'   these are not updated by the shrinkage step.}
#'
#' \item{upper}{Upper limits of estimated HPD intervals. Note that
#'   these are not updated by the shrinkage step.}
#'
#' \item{z}{z-scores for posterior mean LFC estimates.}
#'
#' \item{lpval}{-log10 two-tailed p-values obtained from the
#'   z-scores. When \code{shrink.method = "ash"}, this is \code{NA}, and
#'   the s-values are returned instead (see below).}
#'
#' \item{svalue}{s-values returned by
#'   \code{\link[ashr]{ash}}. s-values are analogous to q-values, but
#'   based on the local false sign rate; see Stephens (2016).}
#'
#' \item{lfsr}{When \code{shrink.method = "ash"} only, this output
#'   contains the estimated local false sign rates.}
#'
#' \item{ash}{When \code{shrink.method = "ash"} only, this output
#'   contains the \code{\link[ashr]{ash}} return value (after removing
#'   the \code{"data"}, \code{"result"} and \code{"call"} list
#'   elements).}
#' 
#' \item{F}{Maximum-likelihood estimates of the Poisson model
#'   parameters.}
#'
#' \item{f0}{Maximum-likelihood estimates of the null model
#'    parameters.}
#' 
#' \item{ar}{A vector containing the Metropolis acceptance ratios
#'   from each MCMC run.}
#' 
#' @references
#' Stephens, M. (2016). False discovery rates: a new deal.
#' \emph{Biostatistics} \bold{18}(2), kxw041.
#' \doi{10.1093/biostatistics/kxw041}
#'
#' Zhu, A., Ibrahim, J. G. and Love, M. I. (2019). Heavy-tailed prior
#' distributions for sequence count data: removing the noise and
#' preserving large differences. \emph{Bioinformatics} \bold{35}(12),
#' 2084–2092.
#'
#' @examples
#' # Perform a differential expression (DE) analysis using the previously
#' # fitted multinomial topic model. Note that the de_analysis call could
#' # take several minutes to complete.
#' \donttest{
#' set.seed(1)
#' data(pbmc_facs)
#' de <- de_analysis(pbmc_facs$fit,pbmc_facs$counts)
#'
#' # Compile the DE analysis results for topic 4 into a table, and
#' # rank the genes by the posterior mean log-fold change, limiting to
#' # DE genes identified with low lfsr ("local false sign rate").
#' k <- 4
#' dat <- data.frame(postmean = de$postmean[,k],
#'                   z        = de$z[,k],
#'                   lfsr     = de$lfsr[,k])
#' rownames(dat) <- with(pbmc_facs$genes,paste(symbol,ensembl,sep = "_"))
#' dat <- subset(dat,lfsr < 0.01)
#' dat <- dat[order(dat$postmean,decreasing = TRUE),]
#'
#' # Genes at the top of this ranking are genes that are much more
#' # highly expressed in the topic compared to other topics.
#' head(dat,n = 10)
#'
#' # The genes at the bottom of the ranking are genes that are much less
#' # expressed in the topic.
#' tail(dat,n = 10)
#'
#' # Create a volcano plot from the DE results for topic 4.
#' volcano_plot(de,k = k,ymax = 50,labels = pbmc_facs$genes$symbol)
#' }
#' 
#' @importFrom utils modifyList
#' @importFrom Matrix rowSums
#' @importFrom Matrix colSums
#' @importFrom stats rnorm
#' @importFrom stats runif
#' @importFrom ashr ash
#'
#' @export
#' 
de_analysis <- function (fit, X, s = rowSums(X), pseudocount = 0.01,
                         fit.method = c("scd","em","mu","ccd","glm"),
                         shrink.method = c("ash","none"), lfc.stat = "le",
                         control = list(), verbose = TRUE, ...) {

  # CHECK AND PROCESS INPUTS
  # ------------------------
  # Check and process input argument "fit".
  if (is.matrix(fit)) {
    L <- fit
    if (any((rowSums(L) - 1) > 1e-15))
      warning("\"fit\" is a matrix but may not be topic proportions matrix; ",
              "rowSums(X) should be a vector of all ones")
    m <- ncol(X)
    k <- ncol(fit)
    F <- matrix(1,m,k)
    rownames(F) <- colnames(X)
    colnames(F) <- colnames(L)
    fit <- init_poisson_nmf(X,F = F,L = L)
    fit <- poisson2multinom(fit)
  }
  if (!(inherits(fit,"poisson_nmf_fit") |
        inherits(fit,"multinom_topic_model_fit")))
    stop("Input \"fit\" should be an object of class \"poisson_nmf_fit\" or ",
         "\"multinom_topic_model_fit\", or a matrix of topic proportions")
  if (inherits(fit,"poisson_nmf_fit"))
    fit <- poisson2multinom(fit)
  
  # Check and process input argument "X".
  if (!((is.numeric(X) & is.matrix(X)) | is.sparse.matrix(X)))
    stop("Input argument \"X\" should be a numeric matrix (a \"matrix\" or ",
         "a \"dgCMatrix\")")
  verify.fit.and.count.matrix(X,fit)
  if (is.matrix(X) & is.integer(X))
    storage.mode(X) <- "double"

  # Get the number of rows (n) and columns (m) in the counts matrix, and
  # the number of groups or topics (k).
  n <- nrow(X)
  m <- ncol(X)
  k <- ncol(fit$F)

  # Check input argument "s".
  verify.positive.vector(s)
  if (length(s) != n)
    stop("Input argument \"s\" should be a vector of positive numbers, ",
         "in which length(s) = nrow(X)")
  
  # Check input argument "pseudocount".
  if (any(pseudocount <= 0))
    stop("Input argument \"pseudocount\" should be a positive number")

  # Process input arguments "fit.method" and "shrink.method".
  fit.method <- match.arg(fit.method)
  shrink.method <- match.arg(shrink.method)

  # Check and process input argument "control".
  control <- modifyList(de_analysis_control_default(),control,keep.null = TRUE)
  if (control$nc > 1 & .Platform$OS.type == "windows")
    stop("Multithreading is not available on Windows; try again with ",
         "control$nc = 1")
  
  # Check and process input argument "lfc.stat".
  if (!(all(lfc.stat == "vsnull") | all(lfc.stat == "le"))) {
    if (!(any(lfc.stat == 1:k) | any(lfc.stat == colnames(fit$F))))
      stop("Input argument \"lfc.stat\" should be either \"vsnull\", \"le\" ",
           "a number between 1 and k, where k is the number of topics, or ",
           "a name of a topic (column of fit$F)")
    if (is.character(lfc.stat))
      lfc.stat <- match(lfc.stat,colnames(fit$F))
  }
    
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
  L <- pmax(L,control$minval)
  L <- pmin(L,1 - control$minval)

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
  nc <- initialize.multithreading(control$nc,verbose)
  F <- fit_poisson_models(X,L,fit.method,control$eps,control$numiter,
                          control$tol,control$nc)
  F <- pmax(F,control$minval)
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
  ns <- control$ns
  D <- matrix(rnorm(ns*k),ns,k)
  U <- matrix(runif(ns*k),ns,k)
  M <- matrix(sample(k,ns*k,replace = TRUE),ns,k) - 1
  if (nc == 1)
    out <- compute_lfc_stats(X,F,L,f0,D,U,M,lfc.stat,control$conf.level,
                             control$rw,control$eps,verbose)
  else {
    out <- compute_lfc_stats_multicore(X,F,L,f0,D,U,M,lfc.stat,
                                       control$conf.level,control$rw,
                                       control$eps,control$nc,control$nsplit,
                                       verbose)
  }
  if (any(out$ar == 0))
    warning("One or more MCMC simulations yielded acceptance rates of zero; ",
            "consider increasing the number of Monte Carlo samples ",
            "(control$ns) or modifying the noise level of the random-walk ",
            "proposal distribution (control$rw) to improve the acceptance ",
            "rates")

  # STABILIZE ESTIMATES USING ADAPTIVE SHRINKAGE
  # --------------------------------------------
  # If requested, use adaptive shrinkage to stabilize the log-fold
  # change estimates. Here we need to carefully edge cases such as
  # se's of zero.
  if (shrink.method == "ash") {
    if (verbose)
      cat("Stabilizing posterior log-fold change estimates using adaptive",
          "shrinkage.\n")
    se <- with(out,postmean/z)
    se[out$z == 0] <- as.numeric(NA)
    se[out$postmean == 0] <- 0
    res          <- shrink_estimates(out$postmean,se,...)
    out$postmean <- res$b
    out$z        <- res$z
    out$lfsr     <- res$lfsr
    out$lpval    <- as.numeric(NA)
    out$svalue   <- res$svalue
    out$ash      <- res$ash
    dimnames(out$lfsr)   <- dimnames(F)
    dimnames(out$svalue) <- dimnames(F)
  } else {

    # Compute the -log10 two-tailed p-values computed from the z-scores.
    out$lpval  <- -lpfromz(out$z)
    out$svalue <- as.numeric(NA)
    out$lfsr   <- as.numeric(NA)
  }

  # Return the Poisson model MLEs (F), the log-fold change statistics
  # (est, postmean, lower, upper, z, lpval) and local false sign
  # rates (lfsr), and the relative rates under the "null" model (f0).
  out$F  <- F
  out$f0 <- f0
  class(out) <- c("topic_model_de_analysis","list")
  return(out)
}

# Perform adaptive shrinkage on the matrix of effect estimates b and
# their standard errors, then output the revised effect estimates (b),
# standard errors (se) and z-scores. All effects i in which either
# b[i] or se[i] is missing (NA) are not revised.
shrink_estimates <- function (b, se, ...) {
  
  # Set up the z-scores output.
  z <- b
  z[is.na(b) | is.na(se)] <- as.numeric(NA)
  
  # Run adaptive shrinkage.
  i   <- which(!(is.na(b) | is.na(se)))
  out <- ash(b[i],se[i],mixcompdist = "normal",method = "shrink",...)
  b1  <- out$result$PosteriorMean
  se1 <- out$result$PosteriorSD
  
  # Extract the posterior estimates and their standard errors. 
  b[i]  <- b1
  se[i] <- se1
  z[i]  <- b[i]/se[i]
  z[i[b1 == 0]]  <- 0
  z[i[se1 == 0]] <- as.numeric(NA)

  # Extract the lfsr estimates.
  m         <- nrow(b)
  k         <- ncol(b)
  lfsr      <- matrix(as.numeric(NA),m,k)
  svalue    <- matrix(as.numeric(NA),m,k)
  lfsr[i]   <- out$result$lfsr
  svalue[i] <- out$result$svalue

  # Output the revised estimates (b), the standard errors (se), the
  # z-scores (z), the local false sign rates (lfsr), the s-values
  # (svalue) and the raw ash output (ash).
  out[c("data","result","call")] <- NULL
  return(list(b = b,se = se,z = z,lfsr = lfsr,svalue = svalue,ash = out))
}

#' @rdname de_analysis
#'
#' @export
#' 
de_analysis_control_default <- function()
  list(numiter    = 20,
       minval     = 1e-10,
       tol        = 1e-8,
       conf.level = 0.68,
       ns         = 1000,
       rw         = 0.3,
       eps        = 1e-15,
       nc         = 1,
       nsplit     = 100)

# Select genes based on a de_analysis result. Input "de" is an object
# of class "topic_model_de_analysis"; k is a topic, and "subset" is a
# function that returns TRUE or FALSE depending on whether to select
# the gene based on the DE analysis result. Here, "rank" and
# "quantile" should be interpreted as the rank and quantile of the
# gene based on the posterior mean estimate of the log-fold change, so
# that largest positive log-fold change has rank 1 and quantile 1.
# The output is a vector of gene indices (e.g., corresponding to
# rows of de$postmean).
select_de_genes <-
  function (de, k,
            subset = function (postmean, lpval, lfsr, rank, quantile)
              lfsr < 0.05) {

  # Compile data used to select genes.
  n <- nrow(de$postmean)
  dat <- data.frame(postmean = de$postmean[,k],
                    lpval    = de$lpval[,k],
                    lfsr     = as.numeric(NA),
                    rank     = rank(-de$postmean[,k],na.last = TRUE,
                                    ties.method = "average"),
                    quantile = rank(de$postmean[,k],na.last = FALSE,
                                    ties.method = "average")/n)
  if (!all(is.na(de$lfsr)))
    dat$lfsr <- de$lfsr[,k]
  
  # Select the genes according to the provided subset function.
  x <- rep(FALSE,n)
  names(x) <- rownames(de$postmean)
  for (i in 1:n)
    x[i] <- all(subset(dat[i,"postmean"],dat[i,"lpval"],dat[i,"lfsr"],
                       dat[i,"rank"],dat[i,"quantile"]))
  return(which(x))
}
