#' @rdname embeddings_from_topics
#'
#' @title Low-dimensional Embeddings from Poisson NMF or Multinomial Topic Model
#' @description Lightweight interface for rapidly producing
#'   low-dimensional embeddings from matrix factorizations or
#'   multinomial topic models. The defaults used are more suitable for
#'   producing embeddings from matrix factorizations or topic models.
#' 
#' @details Note that since \code{tsne_from_topics} and
#'   \code{umap_from_topics} use nonlinear transformations of the data,
#'   distances between points are generally less interpretable than a
#'   linear transformation obtained by, say, PCA.
#'
#' @param fit An object of class \dQuote{poisson_nmf_fit} or
#'   \dQuote{multinom_topic_model_fit}.
#'
#' @param dims The number of dimensions in the embedding. In
#'   \code{tsne_from_topics}, this is passed as argument \dQuote{dims}
#'   to \code{\link[Rtsne]{Rtsne}}. In \code{umap_from_topics}, this is
#'   passed as argument \dQuote{n_components} to
#'   \code{\link[uwot]{umap}}.
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
#'   \code{\link[stats]{prcomp}}, \code{\link[Rtsne]{Rtsne}} or
#'   \code{\link[uwot]{umap}}.
#'
#' @return An n x d matrix containing the embedding, where n is the
#'   number of rows of \code{fit$L}, and \code{d = dims}.
#' 
#' @seealso \code{\link{pca_plot}}, \code{\link{tsne_plot}},
#'   \code{\link{umap_plot}}, \code{\link[stats]{prcomp}},
#'   \code{\link[Rtsne]{Rtsne}}, \code{\link[uwot]{umap}}
#'
#' @references
#' Kobak, D. and Berens, P. (2019). The art of using t-SNE for
#' single-cell transcriptomics. \emph{Nature Communications} \bold{10},
#' 5416. \url{https://doi.org/10.1038/s41467-019-13056-x}
#' 
#' @examples
#' library(ggplot2)
#' library(cowplot)
#' set.seed(1)
#' data(pbmc_facs)
#' 
#' # Get the Poisson NMF and multinomial topic model fit to the PBMC data.
#' fit1 <- multinom2poisson(pbmc_facs$fit)
#' fit2 <- pbmc_facs$fit
#' fit2 <- poisson2multinom(fit1)
#' 
#' # Compute the first two PCs of the loadings matrix (for the topic
#' # model, fit2, the loadings are the topic proportions).
#' Y1 <- pca_from_topics(fit1)
#' Y2 <- pca_from_topics(fit2)
#' subpop <- pbmc_facs$samples$subpop
#' quickplot(Y1[,1],Y1[,2],color = subpop) + theme_cowplot()
#' quickplot(Y2[,1],Y2[,2],color = subpop) + theme_cowplot()
#'
#' # Compute a 2-d embedding of the loadings using t-SNE.
#' Y1 <- tsne_from_topics(fit1)
#' Y2 <- tsne_from_topics(fit2)
#' quickplot(Y1[,1],Y1[,2],color = subpop) + theme_cowplot()
#' quickplot(Y2[,1],Y2[,2],color = subpop) + theme_cowplot()
#'
#' # Compute a 2-d embedding of the loadings using UMAP.
#' Y1 <- umap_from_topics(fit1)
#' Y2 <- umap_from_topics(fit2)
#' quickplot(Y1[,1],Y1[,2],color = subpop) + theme_cowplot()
#' quickplot(Y2[,1],Y2[,2],color = subpop) + theme_cowplot()
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
#' @param pca Whether to perform a PCA processing step in t-SNE or
#'   UMAP; passed as argument \dQuote{pca} to \code{\link[Rtsne]{Rtsne}}
#'   or \code{\link[uwot]{umap}}.
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
#'   \code{\link[Rtsne]{Rtsne}} or \code{\link[uwot]{umap}}.
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
  Y <- out$Y
  rownames(Y) <- rownames(fit$L)
  colnames(Y) <- paste0("tsne",1:dims)
  return(Y)
}

#' @rdname embeddings_from_topics
#'
#' @param n_neighbors Number of nearest neighbours in manifold
#'   approximation; passed as argument \dQuote{n_neighbors} to
#'   \code{\link[uwot]{umap}}.
#'
#' @param metric Distance matrix used to find nearest neighbors;
#'   passed as argument \dQuote{metric} to
#'   \code{\link[uwot]{umap}}.
#'
#' @param scale Scaling to apply to \code{fit$L}; passed as argument
#'   \dQuote{scale} to \code{\link[uwot]{umap}}.
#'
#' @importFrom uwot umap
#' 
#' @export
#' 
umap_from_topics <- function (fit, dims = 2, n_neighbors = 30,
                              metric = "euclidean", scale = "none",
                              pca = NULL, verbose = TRUE, ...) {

  # Check input argument "fit".
  if (!(inherits(fit,"poisson_nmf_fit") |
        inherits(fit,"multinom_topic_model_fit")))
    stop("Input \"fit\" should be an object of class \"poisson_nmf_fit\" or ",
         "\"multinom_topic_model_fit\"")
  n <- nrow(fit$L)
  if (n < 4)
    stop("umap cannot run on a matrix with fewer than 4 rows")

  # Adjust n_neighbors if it is too large for the number of samples.
  n0 <- n - 1
  if (n_neighbors > n0) {
    message(sprintf(paste("n_neighbors automatically changed to %d because",
                          "the original setting of %d was too large for the",
                          "number of samples (%d)"),n0,n_neighbors,n))
    n_neighbors <- n0
  }
  
  # Compute the UMAP embedding.
  Y <- umap(fit$L,n_neighbors = n_neighbors,n_components = dims,pca = pca,
            metric = metric,scale = scale,verbose = verbose,...)
  rownames(Y) <- rownames(fit$L)
  colnames(Y) <- paste0("umap",1:dims)
  return(Y)
}
