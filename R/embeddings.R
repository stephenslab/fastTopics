#' @rdname embeddings_from_topics
#'
#' @title Low-dimensional Embeddings from Poisson NMF or Multinomial Topic Model
#' @description Lightweight interface for rapidly producing
#'   low-dimensional embeddings from matrix factorizations or
#'   multinomial topic models. The defaults used are more suitable for
#'   producing embeddings from matrix factorizations or topic models.
#' 
#' @details Note that since \code{tsne_plot} uses a \emph{nonlinear}
#'   transformation of the data, distances between points are less
#'   interpretable than a linear transformation obtained by, say, PCA.
#'
#' @param fit An object of class \dQuote{poisson_nmf_fit} or
#'   \dQuote{multinom_topic_model_fit}.
#'
#' @param dims The number of dimensions in the embedding. In
#'   \code{tsne_from_topics}, this is passed as argument \dQuote{dims}
#'   to \code{\link[Rtsne]{Rtsne}}.
#'
#' @param center A logical value indicating whether columns of
#'   \code{fit$L} should be zero-centered before performing PCA; passed
#'   as argument \dQuote{center} to \code{\link[stats]{prcomp}}.
#'
#' @param scale. A logical value indicating whether columns of
#'   \code{fit$L} should be scaled to have unit variance prior to
#'   performing PCA; passed as argument \dQuote{scale.} to
#'   \code{\link[stats]{prcomp}}.
#' 
#' @param \dots Additional arguments passed to
#'   \code{\link[stats]{prcomp}} or \code{\link[Rtsne]{Rtsne}}.
#'
#' @return An n x d matrix containing the embedding, where n is the
#'   number of rows of \code{fit$L}, and \code{d = dims}.
#' 
#' @seealso \code{\link{pca_plot}}, \code{\link{tsne_plot}},
#'   \code{\link{umap_plot}}, \code{\link[stats]{prcomp}},
#'   \code{\link[Rtsne]{Rtsne}}, \code{\link[uwot]{umap}}
#'
#' @examples
#' # Add examples here.
#' 
#' @importFrom stats prcomp
#' 
#' @export
#' 
pca_from_topics <- function (fit, dims = 2, center = TRUE,
                             scale. = FALSE, ...) {
  if (!(inherits(fit,"poisson_nmf_fit") |
        inherits(fit,"multinom_topic_model_fit")))
    stop("Input \"fit\" should be an object of class \"poisson_nmf_fit\" or ",
         "\"multinom_topic_model_fit\"")
  Y <- prcomp(fit$L,retx = TRUE,center = center,scale. = scale.,...)$x
  return(Y[,1:dims])
}

#' @rdname embeddings_from_topics
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
#' @references
#' Kobak, D. and Berens, P. (2019). The art of using t-SNE for
#' single-cell transcriptomics. \emph{Nature Communications} \bold{10},
#' 5416. \url{https://doi.org/10.1038/s41467-019-13056-x}
#' 
#' @importFrom Rtsne Rtsne
#' 
#' @export
#' 
tsne_from_topics <- function (fit, dims = 2, pca = FALSE, normalize = FALSE,
                              perplexity = 100, theta = 0.1, max_iter = 1000,
                              eta = 200, check_duplicates = FALSE,
                              verbose = TRUE, ...) {
    
  # Check input argument "fit".
  if (!(inherits(fit,"poisson_nmf_fit") |
        inherits(fit,"multinom_topic_model_fit")))
    stop("Input \"fit\" should be an object of class \"poisson_nmf_fit\" or ",
         "\"multinom_topic_model_fit\"")
  n <- nrow(fit$L)
  if (n < 4)
    stop("Rtsne cannot run on a matrix with fewer than 4 rows")
      
  # Adjust the perplexity if it is too large for the number of samples.
  p0 <- pmax(1,floor((n - 1)/3) - 1)
  if (perplexity > p0) {
    message(sprintf(paste("Perplexity automatically changed to %d because",
                          "the original setting of %d was too large for the",
                          "number of samples (%d)"),p0,perplexity,n))
    perplexity <- p0
  }
  
  # Compute the t-SNE embedding.
  out <- Rtsne(fit$L,dims,pca = pca,normalize = normalize,eta = eta,
               perplexity = perplexity,theta = theta,max_iter = max_iter,
               check_duplicates = check_duplicates,verbose = verbose,...)

  # Return the t-SNE embedding.
  Y <- out$Y
  rownames(Y) <- rownames(fit$L)
  colnames(Y) <- paste0("d",1:dims)
  return(Y)
}

#' @rdname embeddings_from_topics
#'
#' @importFrom uwot umap
#' 
#' @export
#' 
umap_from_topics <- function (fit, dims, ...) {

}
