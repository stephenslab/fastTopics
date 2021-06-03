#' @rdname embeddings_from_topics
#'
#' @title Low-dimensional Embeddings from Poisson NMF or Multinomial Topic Model
#'
#' @description Describe pca_from_topics here.
#'
#' @details Add details here.
#'
#' @param fit Describe input argument "fit" here.
#'
#' @param dims Describe input argument "dims" here.
#'
#' @param center Describe input argument "center" here.
#'
#' @param scale. Describe input argument "scale." here.
#' 
#' @param \dots Describe additional inputs here.
#' 
#' @seealso \code{\link{pca_plot}}, \code{\link{tsne_plot}},
#'   \code{\link{umap_plot}}
#'
#' @importFrom stats prcomp
#' 
#' @export
#' 
pca_from_topics <- function (fit, dims = 2, center = TRUE, scale. = FALSE,
                             ...) {
  if (!(inherits(fit,"poisson_nmf_fit") |
        inherits(fit,"multinom_topic_model_fit")))
    stop("Input \"fit\" should be an object of class \"poisson_nmf_fit\" or ",
         "\"multinom_topic_model_fit\"")
  Y <- prcomp(fit$L,retx = TRUE,center = center,scale. = scale.,...)$x
  return(Y[,1:dims])
}

#' @rdname embeddings_from_topics
#' 
umap_from_topics <- function (fit, dims, ...) {

}

#'
#' #' @title t-SNE from Poisson NMF or Multinomial Topic Model
#'
#' @description Computes a low-dimensional nonlinear embededding of
#'   the data from the estimated loadings or mixture proportions using
#'   the t-SNE nonlinear dimensionality reduction method.
#'
#' @details This is a lightweight interface for rapidly producing
#' t-SNE embeddings from matrix factorizations or multinomial topic
#' models; in particular, \code{tsne_from_topics} replaces the t-SNE
#' defaults with settings that are more suitable for visualizing the
#' structure of a matrix factorization or topic model (e.g., the PCA
#' step in \code{Rtsne} is activated by default, but disabled in
#' \code{tsne_from_topics}). See Kobak and Berens (2019) for guidance
#' on choosing t-SNE settings such as the "perplexity" and learning
#' rate (\code{eta}).
#'
#' Note that since \code{tsne_plot} uses a \emph{nonlinear}
#' transformation of the data, distances between points are less
#' interpretable than a linear transformation visualized using
#' \code{pca_plot} for example.
#'
#' @param fit An object of class \dQuote{poisson_nmf_fit} or
#' \dQuote{multinom_topic_model_fit}.
#'
#' @param dims The number of dimensions in the t-SNE embedding; passed
#'   as argument \dQuote{dims} to \code{\link[Rtsne]{Rtsne}}.
#'
#' @param n The maximum number of rows in the loadings matrix
#'   \code{fit$L} to use; when the loadings matrix has more than
#'   \code{n} rows, the t-SNE embedding is computed on a random
#'   selection of \code{n} rows. An upper limit on the number of rows is
#'   used because the runtime of \code{\link[Rtsne]{Rtsne}} increases
#'   rapidly with the number of rows in the input matrix.
#'
#' @param scaling A numeric vector of length equal to the number of
#'   topics specifying a scaling of the columns of \code{fit$L}; this
#'   re-scaling is performed prior to running t-SNE. The vector should
#'   contain non-negative numbers only. A larger value will increase the
#'   importance, or \dQuote{weight}, of the respective topic in
#'   computing the embedding. When \code{scaling} is \code{NULL}, no
#'   re-scaling is performed. Note that this scaling will have no effect
#'   if \code{normalize = TRUE}.
#' 
#' @param pca Whether to perform a PCA processing step in t-SNE;
#'   passed as argument \dQuote{pca} to \code{\link[Rtsne]{Rtsne}}.
#'
#' @param normalize Whether to normalize the data prior to running
#'   t-SNE; passed as argument \dQuote{normalize} to
#'   \code{\link[Rtsne]{Rtsne}}.
#'
#' @param perplexity t-SNE perplexity parameter, passed as argument
#'   \dQuote{perplexity} to \code{\link[Rtsne]{Rtsne}}. The perplexity
#'   is automatically revised if it is too large; see
#'   \code{\link[Rtsne]{Rtsne}} for more information. 
#'
#' @param theta t-SNE speed/accuracy trade-off parameter; passed as
#'   argument \dQuote{theta} to \code{\link[Rtsne]{Rtsne}}.
#' 
#' @param max_iter Maximum number of t-SNE iterations; passed as
#'   argument \dQuote{max_iter} to \code{\link[Rtsne]{Rtsne}}.
#'
#' @param eta t-SNE learning rate parameter; passed as argument
#'   \dQuote{eta} to \code{\link[Rtsne]{Rtsne}}.
#'
#' @param check_duplicates When \code{check_duplicates = TRUE}, checks
#'   whether there are duplicate rows in \code{fit$L}; passed as argument
#'   \dQuote{check_duplicates} to \code{\link[Rtsne]{Rtsne}}.
#' 
#' @param verbose If \code{verbose = TRUE}, progress updates are
#'   printed; passed as argument \dQuote{verbose} to
#'   \code{\link[Rtsne]{Rtsne}}.
#' 
#' @param \dots Additional arguments passed to \code{\link[Rtsne]{Rtsne}}.
#'
#' @return A list with two list elements: \code{Y}, an n x d matrix
#'   containing the embedding \code{Y} returned by
#'   \code{\link[Rtsne]{Rtsne}}, where n is the number of rows of the
#'   loadings matrix, and \code{d = dims}; \code{rows}, the rows of the
#'   loadings matrix included in the t-SNE embedding.
#'
#' @references
#'
#' Kobak, D. and Berens, P. (2019). The art of using t-SNE for
#' single-cell transcriptomics. \emph{Nature Communications} \bold{10},
#' 5416. \url{https://doi.org/10.1038/s41467-019-13056-x}
#' 
#' @seealso \code{\link[Rtsne]{Rtsne}}
#' 
#' @importFrom Rtsne Rtsne
#' 
#' @export
#' 
tsne_from_topics <- function (fit, dims = 2, n = 5000, scaling = NULL,
                              pca = FALSE, normalize = FALSE,
                              perplexity = 100, theta = 0.1, max_iter = 1000,
                              eta = 200, check_duplicates = FALSE,
                              verbose = TRUE, ...) {
    
  # Check input argument "fit".
  if (!(inherits(fit,"poisson_nmf_fit") |
        inherits(fit,"multinom_topic_model_fit")))
    stop("Input \"fit\" should be an object of class \"poisson_nmf_fit\" or ",
         "\"multinom_topic_model_fit\"")

  # Randomly subsample, if necessary.
  n0 <- nrow(fit$L)
  if (n < n0)
    rows <- sample(n0,n)
  else
    rows <- 1:n0

  # Prepare the data for t-SNE. If requested, re-scale the columns of
  # the loadings matrix.
  L <- fit$L[rows,]
  if (!is.null(scaling))
    L <- scale.cols(L,scaling)

  # Adjust the perplexity if it is too large for the number of samples.
  n  <- nrow(L)
  p0 <- floor((n - 1)/3) - 1
  if (perplexity > p0) {
    message(sprintf(paste("Perplexity automatically changed to %d because",
                          "original setting of %d was too large for the",
                          "number of samples (%d)"),p0,perplexity,n))
    perplexity <- p0
  }
  
  # Compute the t-SNE embedding.
  out <- Rtsne(L,dims,pca = pca,normalize = normalize,perplexity = perplexity,
               theta = theta,max_iter = max_iter,eta = eta,
               check_duplicates = check_duplicates,
               verbose = verbose,...)

  # Return the t-SNE embedding stored as an n x dims matrix (Y), and
  # the rows of L included in the embedding (rows).
  Y           <- out$Y
  rownames(Y) <- rownames(L)
  colnames(Y) <- paste0("d",1:dims)
  return(list(Y = Y,rows = rows))
}

