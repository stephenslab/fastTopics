#' @title PCA Plot
#'
#' @description Visualize the structure of the Poisson NMF loadings or
#'   the multinomial topic model mixture proportions by projection onto
#'   principal components (PCs). \code{plot_hexbin_plot} is most useful
#'   for visualizing the PCs of a data set with thousands of samples, or
#'   more.
#'
#' @details This is a lightweight interface for rapidly producing PCA
#' plots from a Poisson non-negative matrix factorization or
#' multinomial topic model. The plotting is implemented using
#' ggplot2. For more control over the plot's appearance, the plot can
#' be customized by modifying the \code{ggplot_call} and
#' \code{plot_grid_call} arguments. If not provided, the principal
#' components are computed using \code{\link[stats]{prcomp}}.
#'
#' Since principal components are a \emph{linear} projection
#' (\dQuote{rotation}) of the data, distances between points can bee
#' more interpretable than nonlinear embedding methods such as t-SNE
#' (as visualized using \code{tsne_plot}) and UMAP.
#'
#' @param fit An object of class \dQuote{poisson_nmf_fit} or
#'   \dQuote{multinom_topic_model_fit}.
#'
#' @param k The dimensions or topics selected by number or name. When
#'   \code{fill = "loading"}, one plot is created per selected dimension
#'   or topic; when \code{fill = "loading"} and \code{k} is not
#'   specified, all dimensions or topics are plotted.
#'
#' @param out.pca A list containing the result of a principal
#'   components analysis, typically the result of calling
#'   \code{\link[stats]{prcomp}} or a comparable function such as
#'   \code{rpca} from the rsvd package. It should be a list object
#'   containing, at a minimum, a list matrix or data frame
#'   \code{out.pca$x} storing the rotated data. If not provided,
#'   principal components will be computed automatically by calling
#'   \code{\link[stats]{prcomp}}.
#'
#' @param pcs The two principal components (PCs) to be plotted,
#'   specified by name or number.
#'
#' @param fill The quantity to map onto the fill colour of the points
#'   in the PCA plot. Set \code{fill = "loading"} to vary the fill
#'   colour according to the loadings (or mixture proportions) of the
#'   select topic, or topics. Alternatively, \code{fill} may be set to a
#'   data vector with one entry per row of \code{fit$L}, in which case
#'   these data are mapped to the fill colour of the points. When
#'   \code{fill = "none"}, the fill colour is not varied.
#' 
#' @param ggplot_call The function used to create the plot. Replace
#'   \code{pca_plot_ggplot_call} or \code{pca_hexbin_plot_ggplot_call}
#'   with your own function to customize the appearance of the plot.
#'
#' @param plot_grid_call When multiple topics (\code{k}) are selected,
#'   and \code{fill = "loading"}, this is the function used to arrange the
#'   plots into a grid using \code{\link[cowplot]{plot_grid}}. It should
#'   be a function accepting a single argument, \code{plots}, a list of
#'   \code{ggplot} objects.
#' 
#' @param \dots Additional arguments passed to
#'   \code{\link[stats]{prcomp}}. These additional arguments are only
#'   used if \code{out.pca} is not provided.
#' 
#' @return A \code{ggplot} object.
#'
#' @export
#' 
pca_plot <-
  function (fit, k, out.pca, pcs = 1:2, fill = "loading",
            ggplot_call = pca_plot_ggplot_call,
            plot_grid_call = function (plots) do.call(plot_grid,plots),
            ...) {

  # Check and process input "fit".
  if (!(inherits(fit,"poisson_nmf_fit") |
        inherits(fit,"multinom_topic_model_fit")))
    stop("Input \"fit\" should be an object of class \"poisson_nmf_fit\" or ",
         "\"multinom_topic_model_fit\"")

  # Process input "k". 
  if (missing(k)) {
    if (all(fill == "loading"))
      k <- seq(1,ncol(fit$L))
    else
      k <- 1
  }

  # Process input "fill".
  fill.label <- deparse(substitute(fill))
  if (!(is.numeric(fill) | all(fill == "loading") | all(fill == "none")))
    if (!is.factor(fill))
      fill <- factor(fill)
  
  # If necessary, compute the principal components using prcomp.
  if (missing(out.pca))
    out.pca <- prcomp(fit$L,...)

  if (length(k) > 1) {

    # Create a PCA plot for each selected topic, and combine them
    # using plot_grid. This is done by recursively calling pca_plot.
    m     <- length(k)
    plots <- vector("list",m)
    names(plots) <- k
    for (i in 1:m)
      plots[[i]] <- pca_plot(fit,k[i],out.pca,pcs,fill,ggplot_call,
                             plot_grid_call,...)
    return(plot_grid_call(plots))
  } else {
      
    # Get the data (y) mapped to fill colour.
    y <- factor(rep("none",nrow(fit$L)))
    fill.type  <- "none"
    if (length(fill) == 1) {
      if (fill == "loading") {
        y          <- fit$L[,k]
        fill.type  <- "loading"
        fill.label <- paste("topic",k)
      }
    } else {
      y <- fill
      if (is.numeric(y))
        fill.type <- "numeric"
      else if (is.factor(y))
        fill.type <- "factor"
    }

    # Prepare the data for plotting.
    if (is.numeric(pcs))
      pcs <- colnames(out.pca$x)[pcs]
    dat <- cbind(as.data.frame(out.pca$x[,pcs]),y)
    names(dat) <- c(pcs,"y")
    
    # Create the PCA plot.
    return(ggplot_call(dat,pcs,fill.type,fill.label))
  }
}

#' @rdname pca_plot
#'
#' @param dat A data frame passed as input to
#'   \code{\link[ggplot2]{ggplot}}, containing, at a minimum, the
#'   principal components to be plotted. The data frame passed to
#'   \code{pca_plot_ggplot_call} should also have a \dQuote{y} column,
#'   which is mapped to the fill colour of the points in the plot.
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
pca_plot_ggplot_call <-
  function (dat, pcs, fill.type = c("loading","numeric","factor","none"),
            fill.label, font.size = 9) {
  fill.type <- match.arg(fill.type)
  p <- ggplot(dat,aes_string(x = pcs[1],y = pcs[2],fill = "y")) +
    geom_point(shape = 21,color = "white",stroke = 0.3) +
    labs(x = pcs[1],y = pcs[2],fill = fill.label) +
    theme_cowplot(font.size) +
    theme(plot.title = element_text(size = font.size,face = "plain"))
  if (fill.type == "loading")
    p <- p + scale_fill_gradient2(low = "lightskyblue",mid = "gold",
                                  high = "orangered",
                                  midpoint = mean(range(dat$y)))
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

#' @rdname pca_plot
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
pca_hexbin_plot <- function (fit, out.pca, pcs = 1:2, bins = 40,
                             breaks = c(0,1,10,100,1000,Inf),
                             ggplot_call = pca_hexbin_plot_ggplot_call,
                             ...) {
  # Check and process inputs.
  if (!(inherits(fit,"poisson_nmf_fit") |
        inherits(fit,"multinom_topic_model_fit")))
    stop("Input \"fit\" should be an object of class \"poisson_nmf_fit\" or ",
         "\"multinom_topic_model_fit\"")

  # If necessary, compute the principal components using prcomp.
  if (missing(out.pca))
    out.pca <- prcomp(fit$L,...)

  # Prepare the data for plotting.
  dat <- as.data.frame(out.pca$x)
  if (is.numeric(pcs))
    pcs <- names(dat)[pcs]
  
  # Create the PCA plot.
  return(ggplot_call(dat,pcs,bins,breaks))
}

#' @rdname pca_plot
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

#' @title t-SNE Plot
#'
#' @description Visualize the "structure" of the Poisson NMF loadings
#'   or multinomial topic model mixture proportions by projection onto a
#'   2-d surface. Samples in the projection are colored according to the
#'   their loadings/proportions. By default, t-SNE is used to compute
#'   the 2-d embedding from the loadings or mixture proportions.
#'
#' @details This is a lightweight interface primarily intended to
#' expedite creation of scatterplots for visualizing the loadings or
#' mixture proportions in 2-d; most of the \dQuote{heavy lifting} is
#' done by ggplot2. The 2-d embedding itself is computed by invoking
#' function \code{\link{tsne_from_topics}} (unless the \dQuote{tsne}
#' input is provided). For more control over the plot's appearance,
#' the plot can be customized by modifying the \code{ggplot_call} and
#' \code{plot_grid_call} arguments.
#'
#' An effective 2-d visualization may also necessitate some
#' fine-tunning of the t-SNE settings, such as the
#' \dQuote{perplexity}, or the number of samples included in the
#' plot. The t-SNE settings can be controlled by the additional
#' arguments (\dots) passed to \code{tsne_from_topics}; see
#' \code{\link{tsne_from_topics}} for details. Alternatively, a 2-d
#' embedding may be pre-computed, and passed as argument \code{tsne}
#' to \code{tsne_plot}.
#' 
#' @param fit An object of class \dQuote{poisson_nmf_fit} or
#'   \dQuote{multinom_topic_model_fit}.
#'
#' @param color The data mapped to the color \dQuote{aesthetic} in the
#'   plot. When \code{color = "loading"}, the estimated loadings (stored
#'   the \code{L} matrix) in the Poisson NMF model are plotted; when
#'   \code{color = "mixprop"}, the estimated mixture proportions (which are
#'   recovered from the loadings by calling
#'   \code{\link{poisson2multinom}} are shown. When \code{fit} is a
#'   \dQuote{multinom_topic_model_fit} object, the only available option
#'   is \code{color = "mixprop"}. In most settings, the mixture proportions
#'   are preferred, even if the 2-d embedding is computed from the
#'   loading matrix of the Poisson NMF model.
#' 
#' @param k The topic, or topics, selected by number or name. One plot
#'   is created per selected topic. When not specified, all topics are
#'   plotted.
#'
#' @param tsne A 2-d embedding of the samples (rows of X), or a subset
#'   of the samples, such as an output from
#'   \code{\link{tsne_from_topics}}. It should be a list object with the
#'   same structure as a \code{tsne_from_topics} output; see
#'   \code{\link{tsne_from_topics}} for details. If not provided, a 2-d
#'   t-SNE embedding will be estimated automatically by calling
#'   \code{tsne_from_topics}.
#'
#' @param ggplot_call The function used to create the plot. Replace
#'   \code{tsne_plot_ggplot_call} with your own function to
#'   customize the appearance of the plot.
#'
#' @param plot_grid_call When multiple topics are selected, this is
#'   the function used to arrange the plots into a grid using
#'   \code{\link[cowplot]{plot_grid}}. It should be a function accepting
#'   a single argument, \code{plots}, a list of \code{ggplot} objects.
#' 
#' @param \dots Additional arguments passed to
#'   \code{tsne_from_topics}. These arguments are only used if
#'   \code{tsne} is not provided.
#' 
#' @return A \code{ggplot} object.
#'
#' @importFrom cowplot plot_grid
#' 
#' @export
#' 
tsne_plot <-
  function (fit, color = c("mixprop","loading"), k, tsne,
            ggplot_call = tsne_plot_ggplot_call,
            plot_grid_call = function (plots) do.call(plot_grid,plots), ...) {

  # Check and process inputs.
  color <- match.arg(color)
  if (color == "loading") {
    if (!inherits(fit,"poisson_nmf_fit"))
      stop("For color = \"loading\", input \"fit\" should be an object of ",
           "class \"poisson_nmf_fit\"")
  } else if (color == "mixprop") {
    if (!(inherits(fit,"poisson_nmf_fit") |
          inherits(fit,"multinom_topic_model_fit")))
    stop("Input \"fit\" should be an object of class \"poisson_nmf_fit\" or ",
         "\"multinom_topic_model_fit\"")
  }
  if (missing(k))
    k <- seq(1,ncol(fit$L))
  
  # If necessary, compute the 2-d embedding using t-SNE.
  if (missing(tsne))
    tsne <- tsne_from_topics(fit,2,...)
  
  if (length(k) == 1) {
      
    # Prepare the data for plotting. Note that it is important to
    # subset the rows *after* recovering the parameters of the
    # multinomial topic model.
    if (inherits(fit,"poisson_nmf_fit") & color == "mixprop")
      fit <- poisson2multinom(fit)
    dat <- as.data.frame(cbind(tsne$Y,fit$L[tsne$rows,k]))
    names(dat) <- c("d1","d2","loading")
      
    # Create the t-SNE plot.
    return(ggplot_call(dat,k))
  } else {

    # Create a t-SNE plot for each selected topic, and combine them
    # using plot_grid. This is done by recursively calling tsne_plot.
    m     <- length(k)
    plots <- vector("list",m)
    names(plots) <- k
    for (i in 1:m)
      plots[[i]] <- tsne_plot(fit,color,k[i],tsne,ggplot_call,
                              plot_grid_call,...)
    return(plot_grid_call(plots))
  }
}

#' @rdname tsne_plot
#'
#' @param dat A data frame passed as input to
#'   \code{\link[ggplot2]{ggplot}}, containing, at a minimum, columns
#'   \dQuote{d1}, \dQuote{d2} (the first and second dimensions in the
#'   2-d embedding), and \dQuote{loading}.
#'
#' @param topic.label The name or number of the topic being plotted;
#'   it is only used to determine the plot title.
#' 
#' @param font.size Font size used in plot.
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 scale_fill_gradient2
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_text
#' @importFrom cowplot theme_cowplot
#' 
#' @export
#' 
tsne_plot_ggplot_call <- function (dat, topic.label, font.size = 9)
  ggplot(dat,aes_string(x = "d1",y = "d2",fill = "loading")) +
    geom_point(shape = 21,color = "white",stroke = 0.3) +
    scale_fill_gradient2(low = "lightskyblue",mid = "gold",high = "orangered",
                         midpoint = mean(range(dat$loading))) +
    labs(x = "tsne 1",y = "tsne 2",title = paste("topic",topic.label)) +
    theme_cowplot(font.size) +
    theme(plot.title = element_text(size = font.size,face = "plain"))
