#' @rdname embedding_plots
#'
#' @title PCA, t-SNE and UMAP Plots
#'
#' @description Visualize the structure of the Poisson NMF loadings or
#'   the multinomial topic model topic proportions by projection onto
#'   a 2-d surface. \code{plot_hexbin_plot} is most useful for
#'   visualizing the PCs of a data set with thousands of samples or
#'   more.
#'
#' @details This is a lightweight interface primarily intended to
#' expedite creation of plots for visualizing the loadings or topic
#' proportions; most of the heavy lifting is done by ggplot2. The 2-d
#' embedding itself is computed by invoking
#' \code{\link{pca_from_topics}}, \code{\link{tsne_from_topics}} or
#' \code{\link{umap_from_topics}}. For more control over the plot's
#' appearance, the plot can be customized by modifying the
#' \code{ggplot_call} and \code{plot_grid_call} arguments.
#'
#' An effective 2-d visualization may also require some fine-tunning
#' of the settings, such as the t-SNE \dQuote{perplexity}, or the
#' number of samples included in the plot. The PCA, UMAP, t-SNE
#' settings can be controlled by the additional arguments
#' (\dots). Alternatively, a 2-d embedding may be pre-computed, and
#' passed as argument \code{Y}.
#' 
#' @param fit An object of class \dQuote{poisson_nmf_fit} or
#'   \dQuote{multinom_topic_model_fit}.
#'
#' @param Y The n x 2 matrix containing the 2-d embedding, where n is
#'   the number of rows in \code{fit$L}. If not provided, the embedding
#'   will be computed automatically.
#'
#' @param fill The quantity to map onto the fill colour of the points
#'   in the PCA plot. Set \code{fill = "loading"} to vary the fill
#'   colour according to the loadings (or topic proportions) of the
#'   select topiced or topics. Alternatively, \code{fill} may be set to a
#'   data vector with one entry per row of \code{fit$L}, in which case
#'   the data are mapped to the fill colour of the points. When
#'   \code{fill = "none"}, the fill colour is not varied.
#' 
#' @param k The dimensions or topics selected by number or name. When
#'   \code{fill = "loading"}, one plot is created per selected dimension
#'   or topic; when \code{fill = "loading"} and \code{k} is not
#'   specified, all dimensions or topics are plotted.
#'
#' @param fill.label The label used for the fill colour legend.
#' 
#' @param ggplot_call The function used to create the plot. Replace
#'   \code{embedding_plot_2d_ggplot_call} or \code{pca_hexbin_plot_ggplot_call}
#'   with your own function to customize the appearance of the plot.
#'
#' @param plot_grid_call When \code{fill = "loading"} and multiple
#'   topics (\code{k}) are selected, this is the function used to
#'   arrange the plots into a grid using \code{\link[cowplot]{plot_grid}}.
#'   It should be a function accepting a single argument, \code{plots},
#'   a list of \code{ggplot} objects.
#' 
#' @return A \code{ggplot} object.
#'
#' @seealso \code{\link{pca_from_topics}},
#'   \code{\link{tsne_from_topics}},
#'   \code{\link{umap_from_topics}}
#'
#' @examples
#' set.seed(1)
#' data(pbmc_facs)
#' 
#' # Get the Poisson NMF and multinomial topic models fitted to the
#' # PBMC data.
#' fit1 <- multinom2poisson(pbmc_facs$fit)
#' fit2 <- pbmc_facs$fit
#' fit2 <- poisson2multinom(fit1)
#'
#' # Plot the first two PCs of the loadings matrix (for the
#' # multinomial topic model, "fit2", the loadings are the topic
#' # proportions).
#' subpop <- pbmc_facs$samples$subpop
#' p1 <- pca_plot(fit1,k = 1)
#' p2 <- pca_plot(fit2)
#' p3 <- pca_plot(fit2,fill = "none")
#' p4 <- pca_plot(fit2,pcs = 3:4,fill = "none")
#' p5 <- pca_plot(fit2,fill = fit2$L[,1])
#' p6 <- pca_plot(fit2,fill = subpop)
#'
#' # Plot the loadings using t-SNE.
#' p1 <- tsne_plot(fit1,k = 1)
#' p2 <- tsne_plot(fit2)
#' p3 <- tsne_plot(fit2,fill = subpop)
#'
#' # Plot the loadings using UMAP.
#' p1 <- umap_plot(fit1,k = 1)
#' p2 <- umap_plot(fit2)
#' p3 <- umap_plot(fit2,fill = subpop)
#' 
#' @importFrom cowplot plot_grid
#' 
#' @export
#'
embedding_plot_2d <-
  function (fit, Y, fill = "loading", k, fill.label,
            ggplot_call = embedding_plot_2d_ggplot_call,
            plot_grid_call = function (plots) do.call(plot_grid,plots)) {
    
  # Check input "fit".
  if (!(inherits(fit,"poisson_nmf_fit") |
        inherits(fit,"multinom_topic_model_fit")))
    stop("Input \"fit\" should be an object of class \"poisson_nmf_fit\" or ",
         "\"multinom_topic_model_fit\"")

  # Process input argument Y.
  if (!(is.matrix(Y) & nrow(Y) == nrow(fit$L) & ncol(Y) == 2))
    stop("Input argument \"Y\" should be a matrix with 2 columns and with ",
         "the same number of rows as fit$L")
  if (is.null(colnames(Y)))
    colnames(Y) <- c("d1","d2")

  # Process inputs "fill", "fill.label" and "k".
  if (missing(fill.label))
    fill.label <- deparse(substitute(fill))
  if (!(is.numeric(fill) | all(fill == "loading") | all(fill == "none")))
    if (!is.factor(fill))
      fill <- factor(fill)
  if (missing(k)) {
    if (all(fill == "loading"))
      k <- seq(1,ncol(fit$L))
    else
      k <- 1
  }
  
  if (all(fill == "loading") & length(k) > 1) {

    # Create a plot for each selected dimension (topic) k, and combine them
    # using plot_grid. This is done by recursively calling embedding_plot_2d.
    m     <- length(k)
    plots <- vector("list",m)
    names(plots) <- k
    for (i in 1:m)
      plots[[i]] <- embedding_plot_2d(fit,Y,fill,k[i],fill.label,ggplot_call,
                                      plot_grid_call)
    return(plot_grid_call(plots))
  } else {
      
    # Get the data (y) mapped to the fill colour.
    fill.type <- "none"
    if (length(fill) == 1) {
      if (fill == "loading") {
        fill       <- fit$L[,k]
        fill.type  <- "loading"
        fill.label <- paste("topic",k)
      }
    } else {
      if (is.numeric(fill))
        fill.type <- "numeric"
      else if (is.factor(fill))
        fill.type <- "factor"
    }
  }

  # Create the embedding plot.
  return(embedding_plot_2d_ggplot_call(Y,fill,fill.type,fill.label))
}

#' @rdname embedding_plots
#'
#' @param fill.type The type of variable mapped to fill colour. The
#'   fill colour is not varied when \code{fill.type = "none"}.
#' 
#' @param fill.label The label used for the fill colour legend.
#' 
#' @param font.size Font size used in plot.
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggplot2 scale_fill_gradient2
#' @importFrom ggplot2 scale_fill_gradientn
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_text
#' @importFrom cowplot theme_cowplot
#' 
#' @export
#' 
embedding_plot_2d_ggplot_call <-
  function (Y, fill, fill.type = c("loading","numeric","factor","none"),
            fill.label, font.size = 9) {
  dims <- colnames(Y)
  fill.type <- match.arg(fill.type)
  dat <- cbind(Y,data.frame(fill = fill))
  p <- ggplot(dat,aes_string(x = dims[1],y = dims[2],fill = "fill")) +
    geom_point(shape = 21,color = "white",stroke = 0.3) +
    labs(x = dims[1],y = dims[2],fill = fill.label) +
    theme_cowplot(font.size) +
    theme(plot.title = element_text(size = font.size,face = "plain"))
  if (fill.type == "loading")
    p <- p + scale_fill_gradient2(low = "lightskyblue",mid = "gold",
                                  high = "orangered",
                                  midpoint = mean(range(dat$fill)))
  else if (fill.type == "numeric")
    p <- p + scale_fill_gradientn(na.value = "lightskyblue",
                                  colors = c("skyblue","gold","darkorange",
                                             "magenta"))
  else if (fill.type == "factor")
    p <- p + scale_fill_manual(values = c("#e41a1c","#377eb8","#4daf4a",
                                          "#984ea3","#ff7f00","#ffff33",
                                          "#a65628","#f781bf","#999999"),
                               drop = FALSE)
  else if (fill.type == "none")
    p <- p + scale_fill_manual(values = "black",guide = "none")
  return(p)
}

#' @rdname embedding_plots
#'
#' @param pcs The two principal components (PCs) to be plotted,
#'   specified by name or number.
#'
#' @param n The maximum number of points to plot. If \code{n} is less
#'   than the number of rows of \code{fit$L}, the rows are subsampled at
#'   random. This argument is ignored if \code{Y} is provided.
#' 
#' @param \dots Additional arguments passed to
#'   \code{\link{pca_from_topics}}, \code{\link{tsne_from_topics}} or
#'   \code{\link{umap_from_topics}}. These additional arguments are only
#'   used if \code{Y} is not provided.
#' 
#' @importFrom stats prcomp
#' @importFrom cowplot plot_grid
#' 
#' @export
#' 
pca_plot <-
  function (fit, Y, pcs = 1:2, n = 1e4, fill = "loading", k, fill.label,
            ggplot_call = embedding_plot_2d_ggplot_call,
            plot_grid_call = function (plots) do.call(plot_grid,plots), ...) {
  if (!(inherits(fit,"poisson_nmf_fit") |
        inherits(fit,"multinom_topic_model_fit")))
    stop("Input \"fit\" should be an object of class \"poisson_nmf_fit\" or ",
         "\"multinom_topic_model_fit\"")
  if (missing(fill.label))
    fill.label <- deparse(substitute(fill))
  if (missing(Y)) {
    n0 <- nrow(fit$L)
    if (n < n0) {
      rows <- sample(n0,n)
      fit <- select_loadings(fit,rows)
      if (length(fill) == n0)
        fill <- fill[rows]
    }
    Y <- pca_from_topics(fit,dims = ncol(fit$L),...)
  }
  return(embedding_plot_2d(fit,Y[,pcs],fill,k,fill.label,ggplot_call,
                           plot_grid_call))
}

#' @rdname embedding_plots
#'
#' @importFrom cowplot plot_grid
#' 
#' @export
#' 
tsne_plot <-
  function (fit, Y, n = 2000, fill = "loading", k, fill.label,
            ggplot_call = embedding_plot_2d_ggplot_call,
            plot_grid_call = function (plots) do.call(plot_grid,plots), ...) {
  if (!(inherits(fit,"poisson_nmf_fit") |
        inherits(fit,"multinom_topic_model_fit")))
    stop("Input \"fit\" should be an object of class \"poisson_nmf_fit\" or ",
         "\"multinom_topic_model_fit\"")
  if (missing(fill.label))
    fill.label <- deparse(substitute(fill))
  if (missing(Y)) {
    n0 <- nrow(fit$L)
    if (n < n0) {
      rows <- sample(n0,n)
      fit <- select_loadings(fit,rows)
      if (length(fill) == n0)
        fill <- fill[rows]
    }
    Y <- tsne_from_topics(fit,dims = 2,...)
  }
  return(embedding_plot_2d(fit,Y,fill,k,fill.label,ggplot_call,plot_grid_call))
}

#' @rdname embedding_plots
#'
#' @importFrom cowplot plot_grid
#' 
#' @export
#' 
umap_plot <-
  function (fit, Y, n = 2000, fill = "loading", k, fill.label,
            ggplot_call = embedding_plot_2d_ggplot_call,
            plot_grid_call = function (plots) do.call(plot_grid,plots), ...) {
  if (!(inherits(fit,"poisson_nmf_fit") |
        inherits(fit,"multinom_topic_model_fit")))
    stop("Input \"fit\" should be an object of class \"poisson_nmf_fit\" or ",
         "\"multinom_topic_model_fit\"")
  if (missing(fill.label))
    fill.label <- deparse(substitute(fill))
  if (missing(Y)) {
    n0 <- nrow(fit$L)
    if (n < n0) {
      rows <- sample(n0,n)
      fit <- select_loadings(fit,rows)
      if (length(fill) == n0)
        fill <- fill[rows]
    }
    Y <- umap_from_topics(fit,dims = 2,...)
  }
  return(embedding_plot_2d(fit,Y,fill,k,fill.label,ggplot_call,plot_grid_call))
}

#' @rdname embedding_plots
#'
#' @param bins Number of bins used to create hexagonal 2-d
#'   histogram. Passed as the \dQuote{bins} argument to
#'   \code{\link[ggplot2]{stat_bin_hex}}.
#'
#' @param breaks To produce the hexagonal histogram, the counts are
#'   subdivided into intervals based on \code{breaks}. Passed as the
#'   \dQuote{breaks} argument to \code{\link{cut}}.
#'
#' @export
#' 
pca_hexbin_plot <- function (fit, pcs = 1:2, Y, bins = 40,
                             breaks = c(0,1,10,100,1000,Inf),
                             ggplot_call = pca_hexbin_plot_ggplot_call, ...) {
  # Check and process inputs.
  if (!(inherits(fit,"poisson_nmf_fit") |
        inherits(fit,"multinom_topic_model_fit")))
    stop("Input \"fit\" should be an object of class \"poisson_nmf_fit\" or ",
         "\"multinom_topic_model_fit\"")
  if (missing(out.pca))
    out.pca <- prcomp(fit$L,...)
  dat <- as.data.frame(out.pca$x)
  if (is.numeric(pcs))
    pcs <- names(dat)[pcs]
  return(ggplot_call(dat,pcs,bins,breaks))
}

#' @rdname embedding_plots
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 aes_q
#' @importFrom ggplot2 after_stat
#' @importFrom ggplot2 stat_bin_hex
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_text
#' @importFrom cowplot theme_cowplot
#' 
#' @export
#' 
pca_hexbin_plot_ggplot_call <- function (dat, pcs, bins, breaks, font.size = 9)
  ggplot(dat,aes_string(x = pcs[1],y = pcs[2])) +
    stat_bin_hex(mapping = aes_q(fill = quote(cut(after_stat(count),breaks))),
                                 bins = bins) +
    scale_fill_manual(values = c("gainsboro","lightskyblue","gold","orange",
                                 "magenta")) +
    labs(x = pcs[1],y = pcs[2],fill = "count") +
    theme_cowplot(font_size = font.size) +
    theme(plot.title = element_text(size = font.size,face = "plain"))


