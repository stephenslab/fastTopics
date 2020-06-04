#' @title Plot Progress of Poisson NMF Optimization Method Over Time
#'
#' @description Create a plot showing improvement in one or more
#'   Poisson NMF model fits over time. The horizontal axis shows the
#'   recorded runtime (in s), and the vertical axis shows some quantity
#'   measuring the quality of the fit: the log-likelihood, deviance or
#'   maximum residual of the Karush-Kuhn-Tucker (KKT) first-order
#'   optimality conditions. To better visualize log-likelihoods and
#'   deviances, the log-likelihood and deviance differences are shown on
#'   the logarithmic scale. where the differences are calculated with
#'   respect to the best value achieved over all the fits compared.
#'
#' @details Note that only minimal argument checking is performed.
#' 
#' @param fits An object of class \code{"poisson_nmf_fit"}, or a
#'   non-empty, named list in which each list element is an object of
#'   class \code{"poisson_nmf_fit"}.
#'
#' @param x Choose \code{"timing"} to plot improvement in the solution
#'   over time, or choose \code{"iter"} to plot improvement in the
#'   solution per iteration.
#' 
#' @param y Column of the "progress" data frame used to assess
#'   progress of the Poisson NMF optimization method(s). Should be one
#'   of \code{"loglik"} (log-likelihood), \code{"dev"} (deviance) or
#'   \code{"res"} (maximum residual of KKT conditions).
#'
#' @param add.point.every A positive integer giving the iteration
#'   interval for drawing points on the progress curves. Set to
#'   \code{Inf} to prevent points from being drawn on the plot.
#' 
#' @param colors Colours used to draw progress curves; passed as the
#'   \code{values} input to \code{\link[ggplot2]{scale_color_manual}}.
#'   If fewer colours than "fits" are given, the colours are recycled.
#'
#' @param linetypes Line types used to draw progress curves; passed as
#'   the \code{values} input to \code{\link[ggplot2]{scale_linetype_manual}}.
#'   If fewer line types than "fits" are given, the line types are
#'   recycled.
#'
#' @param linesizes Line sizes used to draw progress curves; passed as
#'   the \code{values} input to \code{\link[ggplot2]{scale_size_manual}}.
#'   If fewer line sizes than "fits" are given, the line sizes are
#'    recycled.
#' 
#' @param shapes Shapes used to draw points at the selected
#'   iterations; passed as the \code{values} input to
#'   \code{\link[ggplot2]{scale_shape_manual}}. If fewer shapes than
#'   "fits" are given, the shapes are recycled.
#' 
#' @param fills Fill colours used to draw points at the selected
#'   iterations; passed as the \code{values} input to
#'   \code{\link[ggplot2]{scale_fill_manual}}. If fewer fill colours
#'   than "fits", are given, the fill colours are recycled.
#' 
#' @param e A small, positive number added to the vertical axis (for
#'   \code{y = "loglik"} and \code{y = "dev"} only) so that the
#'   logarithmic scale does not over-emphasize very small differences.
#'
#' @param theme The ggplot2 "theme".
#'
#' @return A \code{ggplot} object.
#'
#' @seealso \code{\link{fit_poisson_nmf}}
#' 
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggplot2 scale_linetype_manual
#' @importFrom ggplot2 scale_size_manual
#' @importFrom ggplot2 scale_shape_manual
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggplot2 labs
#' @importFrom cowplot theme_cowplot
#' 
#' @export
#' 
plot_progress_poisson_nmf <-
  function (fits, x = c("timing","iter"), y = c("loglik","dev","res"),
            add.point.every = 20,
            colors = c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2",
                       "#D55E00","#CC79A7"),
            linetypes = "solid", linesizes = 0.5, shapes = 19, fills = "white",
            e = 0.01, theme = function() theme_cowplot(font_size = 12)) {

  # CHECK & PROCESS INPUTS
  # ----------------------
  # Check and process input "fits". It should either be an object of
  # class poisson_nmf_fit, or a list of poisson_nmf_fit objects.
  if (inherits(fits,"poisson_nmf_fit")) {
    fit.name    <- deparse(substitute(fits))          
    fits        <- list(fits)
    names(fits) <- fit.name
  } else {
    msg <- paste("Input argument \"fits\" should either be an object of",
                 "class \"poisson_nmf_fit\", or a non-empty, named list in",
                 "which each list element is an object of class",
                 "\"poisson_nmf_fit\"")
    if (!(is.list(fits) & !is.null(names(fits)) & length(fits) > 0))
      stop(msg)
    if (!all(sapply(fits,function (x) inherits(x,"poisson_nmf_fit"))))
      stop(msg)
    if (!all(nchar(names(fits)) > 0))
      stop(msg)
  }

  # Check and process input arguments "x" and "y".
  x <- match.arg(x)
  y <- match.arg(y)

  # Check and process input arguments "colors", "linetypes",
  # "linesizes" and "shapes".
  n <- length(fits)
  if (length(colors) < n)
    colors <- rep(colors,length.out = n)
  if (length(linetypes) < n)
    linetypes <- rep(linetypes,length.out = n)
  if (length(linesizes) < n)
    linesizes <- rep(linesizes,length.out = n)
  if (length(shapes) < n)
    shapes <- rep(shapes,length.out = n)
  if (length(fills) < n)
    fills <- rep(fills,length.out = n)

  # PREPARE DATA FOR PLOT
  # ---------------------
  # Combine the progress information from all the fits into one data
  # frame.
  pdat <- prepare_progress_plot_data(fits,e)

  # CREATE PLOT
  # -----------
  # Create the plot showing the improvement in the log-likelihood,
  # deviance, or maximum KKT residual over time.
  return(create_progress_plot(pdat,x,y,add.point.every,colors,linetypes,
                              linesizes,shapes,fills,theme))
}

# Used by plot_progress_poisson_nmf to create a data frame suitable
# for plotting with ggplot.
prepare_progress_plot_data <- function (fits, e) {
  n     <- length(fits)
  labels <- names(fits)
  for (i in 1:n) {
    y         <- fits[[i]]$progress
    y         <- cbind(data.frame("method" = labels[i]),y)
    y$timing  <- cumsum(y$timing)
    fits[[i]] <- y
  }
  out        <- do.call(rbind,fits)
  out$method <- factor(out$method,labels)
  out$loglik <- max(out$loglik) - out$loglik + e
  out$dev    <- out$dev - min(out$dev) + e
  return(out)
}

# Used by plot_progress_poisson_nmf to create the plot.
create_progress_plot <- function (pdat, x, y, add.point.every, colors,
                                  linetypes, linesizes, shapes, fills, theme) {
  rows <- which(pdat$iter %% add.point.every == 1)
  if (x == "timing")
    xlab <- "runtime (s)"
  else if (x == "iter")
    xlab <- "iteration"
  if (y == "res")
    ylab <- "max KKT residual"
  else if (y == "loglik" | y == "dev")
    ylab <- paste("distance from best",y)
  return(ggplot(pdat,aes_string(x = x,y = y,color = "method",
                                linetype = "method",size = "method")) +
         geom_line(na.rm = TRUE) +
         geom_point(data = pdat[rows,],
                    mapping = aes_string(x = x,y = y,color = "method",
                                         fill = "method",shape = "method"),
                    inherit.aes = FALSE,na.rm = TRUE) +
         scale_y_continuous(trans = "log10") +
         scale_color_manual(values = colors) +
         scale_linetype_manual(values = linetypes) +
         scale_size_manual(values = linesizes) +
         scale_shape_manual(values = shapes) +
         scale_fill_manual(values = fills) +
         labs(x = xlab,y = ylab) +
         theme())
}

#' @rdname loadings_plot
#' 
#' @title Loadings Plot
#'
#' @description Generate one or more bar plots to visualize the
#'   relationship between the loadings, or topic probabilities, and a
#'   selected categorical variable (a factor).
#'
#' @param fit An object of class \dQuote{poisson_nmf_fit} or
#'   \dQuote{multinom_topic_model_fit}.
#'
#' @param x A categorical variable represented as a
#'   \code{\link{factor}}. It should have the same number of elements as
#'   the number of rows in \code{fit$L}.
#'
#' @param k The topic, or topics, selected by number or name. When not
#'   specified, all topics are plotted.
#' 
#' @param ggplot_call The function used to create the plot. Replace
#'   \code{loadings_plot_ggplot_call} with your own function to
#'   customize the appearance of the plot.
#'
#' @param plot_grid_call When multiple topics are selected, this is
#'   the function used to arrange the plots into a grid using
#'   \code{\link[cowplot]{plot_grid}}. It should be a function accepting
#'   a single argument, \code{plots}, a list of \code{ggplot} objects.
#' 
#' @return A \code{ggplot} object.
#'
#' @importFrom cowplot plot_grid
#' 
#' @export
#' 
loadings_plot <-
  function (fit, x, k, ggplot_call = loadings_plot_ggplot_call,
            plot_grid_call = function (plots) do.call(plot_grid,plots)) {

  # Check and process input arguments.
  if (!(inherits(fit,"poisson_nmf_fit") |
        inherits(fit,"multinom_topic_model_fit")))
    stop("Input \"fit\" should be an object of class \"poisson_nmf_fit\" or ",
         "multinom_topic_model_fit")
  n <- nrow(fit$L)
  if (!(is.factor(x) && length(x) == n))
    stop("Input \"x\" should be a factor with as many elements as there ",
         "are rows in the loadings matrix")
  if (missing(k))
    k <- seq(1,ncol(fit$L))
  
  if (length(k) == 1)

    # Create the loadings plot.
    return(ggplot_call(data.frame(x = x,loading = fit$L[,k]),k))
  else {

    # Create a loadings plot for each selected topic, and combine them
    # using plot_grid. This is done by recursively calling loadings_plot.
    m     <- length(k)
    plots <- vector("list",m)
    names(plots) <- k
    for (i in 1:m)
      plots[[i]] <- loadings_plot(fit,x,k[i],ggplot_call)
    return(plot_grid_call(plots))
  }
}

#' @rdname loadings_plot
#'
#' @param dat A data frame passed as input to
#'   \code{\link[ggplot2]{ggplot}}, containing, at a minimum, columns
#'   "x" and "loading".
#'
#' @param topic.label The name or number of the topic being plotted;
#' it is only used to determine the plot title.
#' 
#' @param font.size Font size used in plot.
#' 
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 geom_boxplot
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 element_text
#' @importFrom cowplot theme_cowplot
#' 
#' @export
#' 
loadings_plot_ggplot_call <- function (dat, topic.label, font.size = 10)
  ggplot(dat,aes_string(x = "x",y = "loading")) +
    geom_boxplot(width = 0.25,size = 0.4,outlier.shape = NA) +
    labs(x = "",y = "loading",title = paste("topic",topic.label)) +
    theme_cowplot(font.size) +
    theme(axis.line   = element_blank(),
          axis.text.x = element_text(angle = 45,hjust = 1),
          plot.title  = element_text(size = font.size,face = "plain"))

#' @title t-SNE from Poisson NMF or Multinomial Topic Model
#'
#' @description Computes a low-dimensional embededding of the data from
#'   the estimated loadings, or topic probabilities, using the t-SNE
#'   nonlinear dimensionality reduction method.
#'
#' @param fit An object of class \dQuote{poisson_nmf_fit} or
#'   \dQuote{multinom_topic_model_fit}.
#'
#' @param dims The number of dimensions in the t-SNE embedding; passed
#'   as argument \dQuote{dims} to \code{\link[Rtsne]{Rtsne}}.
#'
#' @param n The maximum number of rows in the loadings matrix
#'   \code{fit$L} to use; when the loadings matrix has more than
#'   \code{n} rows, the t-SNE embedding is computed on a random
#'   selection of \code{n} rows.
#' 
#' @param pca Whether to perform a PCA processing stepe in t-SNE;
#'   passed as argument \dQuote{pca} to \code{\link[Rtsne]{Rtsne}}.
#'
#' @param normalize Whether to normalize the data prior to running
#'   t-SNE; passed as argument \dQuote{normalize} to
#'   \code{\link[Rtsne]{Rtsne}}.
#'
#' @param max_iter Maximum number of t-SNE iterations; passed as
#'   argument \dQuote{max_iter} to \code{\link[Rtsne]{Rtsne}}.
#'
#' @param verbose If \code{verbose = TRUE}, progress updates are
#'   printed; passed as argument \dQuote{verbose} to
#'   \code{\link[Rtsne]{Rtsne}}.
#' 
#' @param \dots Additional arguments passed to
#'   \code{\link[Rtsne]{Rtsne}}.
#'
#' @return A list with two list elements: \code{Y}, an n x d matrix
#'   containing the embedding \code{Y} returned by
#'   \code{\link[Rtsne]{Rtsne}}, where n is the number of rows of the
#'   loadings matrix, and \code{d = dims}; \code{rows}, the rows of the
#'   loadings matrix included in the t-SNE embedding.
#'
#' @seealso \code{\link[Rtsne]{Rtsne}}
#' 
#' @importFrom Rtsne Rtsne
#' 
#' @export
#' 
tsne_from_topics <- function (fit, dims = 2, n = 5000, pca = FALSE,
                              normalize = FALSE, perplexity = 100,
                              theta = 0.1, max_iter = 1000,
                              verbose = TRUE, ...) {
    
  # Check and process input arguments.
  if (!(inherits(fit,"poisson_nmf_fit") |
        inherits(fit,"multinom_topic_model_fit")))
    stop("Input \"fit\" should be an object of class \"poisson_nmf_fit\" or ",
         "multinom_topic_model_fit")

  # Randomly subsample, if necessary.
  L  <- fit$L
  n0 <- nrow(L)
  if (n < n0) {
    rows <- sample(n0,n)
    L    <- L[rows,]
  } else
    rows <- 1:n0

  # Compute the t-SNE embedding.
  out <- Rtsne(L,dims,pca = pca,normalize = normalize,perplexity = perplexity,
               theta = theta,max_iter = max_iter,verbose = verbose,...)

  # Return the t-SNE embedding stored as an n x dims matrix (Y), and
  # the rows of L included in the embedding.
  Y           <- out$Y
  rownames(Y) <- rownames(L)
  colnames(Y) <- paste0("d",1:dims)
  return(list(Y = Y,rows = rows))
}

#' @title t-SNE Plot
#'
#' @description Visualize the "structure" of the Poisson NMF loadings
#'   ("activations") or the topic probabilities by projection onto a 2-d
#'   surface. Samples in the projection are colored according to the their
#'   loadings/probabilities. By default, t-SNE is used to compute the
#'   2-d embedding from the loadings or topic probabilities.
#'
#' @param fit An object of class \dQuote{poisson_nmf_fit} or
#'   \dQuote{multinom_topic_model_fit}.
#'
#' @param color The data mapped to the color \dQuote{aesthetic} in the
#'   plot. When \code{color = "loadings"}, the estimated loadings
#'   (stored the \code{L} matrix) in the Poisson NMF model are plotted;
#'   when \code{color = "probs"}, the estimated topic probabilities
#'   (which can be recovered from the loadings by calling
#'   \code{\link{poisson2multinom}} are shown.
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
#' @param \dots Additional arguments passed to \code{tsne_from_topics}.
#' 
#' @return A \code{ggplot} object.
#'
#' @importFrom cowplot plot_grid
#' 
#' @export
#' 
tsne_plot <-
  function (fit, color = c("probs","loadings"), k, tsne,
            ggplot_call = tsne_plot_ggplot_call,
            plot_grid_call = function (plots) do.call(plot_grid,plots),
            ...) {

  # Check and process inputs.
  color <- match.arg(color)
  if (color == "loadings") {
    if (!inherits(fit,"poisson_nmf_fit"))
      stop("For color = \"loadings\", input \"fit\" should be an object of ",
           "class \"poisson_nmf_fit\"")
  } else if (color == "probs") {
    if (!(inherits(fit,"poisson_nmf_fit") |
          inherits(fit,"multinom_topic_model_fit")))
    stop("Input \"fit\" should be an object of class \"poisson_nmf_fit\" or ",
         "multinom_topic_model_fit")
  }
  if (missing(k))
    k <- seq(1,ncol(fit$L))
  
  # If necessary, compute the t-SNE embedding.
  if (missing(tsne))
    tsne <- tsne_from_topics(fit,2,...)
  
  if (length(k) == 1) {
    if (inherits(fit,"poisson_nmf_fit") & color == "probs")
      fit <- poisson2multinom(fit)
      
    # Prepare the data for plotting.
    fit        <- select(fit,loadings = tsne$rows)
    dat        <- as.data.frame(cbind(tsne$Y,fit$L[,k]))
    names(dat) <- c("d1","d2","loading")
      
    # Create the t-SNE plot.
    return(ggplot_call(dat,k))
  } else {

    # Create a t-SNE plot for each selected topic, and combine them
    # using plot_grid. This is done by recursively calling loadings_plot.
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
#'   "d1", "d2" (the first and second dimensions in the 2-d embedding),
#'   and "loading".
#'
#' @param topic.label The name or number of the topic being plotted;
#' it is only used to determine the plot title.
#' 
#' @param font.size Font size used in plot.
#' 
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 scale_fill_gradient2
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 element_text
#' @importFrom cowplot theme_cowplot
#' 
#' @export
#' 
tsne_plot_ggplot_call <- function (dat, topic.label, font.size = 10)
  ggplot(dat,aes_string(x = "d1",y = "d2",fill = "loading")) +
    geom_point(shape = 21,color = "white",stroke = 0.3) +
    scale_fill_gradient2(low = "lightskyblue",mid = "gold",high = "orangered",
                         midpoint = mean(range(dat$loading))) +
    labs(x = "tsne 1",y = "tsne 2",title = paste("topic",topic.label)) +
    theme_cowplot(font.size) +
    theme(axis.line  = element_blank(),
          plot.title = element_text(size = font.size,face = "plain"))

#' Structure Plot
#'
#' @description Describe the structure plot. Explain why it is called
#' a "structure plot" (add reference(s) below).
#' 
#' @param fit Describe input argument "fit" here.
#'
#' @param topics Describe input argument "topics" here.
#' 
#' @param rows Describe input argument "rows" here.
#'
#' @param colors Describe input argument "colors" here.
#'
#' @param ggplot_call Describe input argument "ggplot_call" here.
#' 
#' @references
#'
#' [Add reference(s) here.]
#' 
#' @export
#'
structure_plot <-
  function (fit, topics, rows,
            colors = c("lightcyan","indianred","gold","dodgerblue",
                       "seagreen","darkmagenta"),
            ggplot_call = structure_plot_ggplot_call) {

  # Check and process inputs.
  if (!(inherits(fit,"poisson_nmf_fit") |
        inherits(fit,"multinom_topic_model_fit")))
    stop("Input \"fit\" should be an object of class \"poisson_nmf_fit\" or ",
         "multinom_topic_model_fit")
  if (inherits(fit,"poisson_nmf_fit"))
    fit <- poisson2multinom(fit)

  # Prepare the data for plotting.
  fit <- select(fit,loadings = rows)
  dat <- compile_structure_plot_data(fit,topics)
      
  # Create the structure plot.
  return(ggplot_call(dat,colors))
}

#' @rdname structure_plot
#'
#' @param dat Describe input argument "dat" here.
#'
#' @param font.size Describe input argument "font.size" here.
#' 
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 geom_col
#' @importFrom ggplot2 scale_x_continuous
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_blank
#' @importFrom cowplot theme_cowplot
#' 
#' @export
#'
structure_plot_ggplot_call <- function (dat, colors, font.size = 10)
  ggplot(dat,aes_string(x = "sample",y = "prob",color = "topic",
                        fill = "topic")) +
    geom_col() +
    scale_x_continuous(breaks = NULL,limits = c(0,max(dat$sample) + 1)) +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors) +
    labs(x = "",y = "topic probability") +
    theme_cowplot(font.size) +
    theme(axis.line = element_blank())

# TO DO: Explain here what this function does, and how to use it.
compile_structure_plot_data <- function (fit, topics) {
  n   <- nrow(fit$L)
  k   <- length(topics)
  dat <- data.frame(sample = rep(1:n,times = k),
                    topic  = rep(topics,each = n),
                    prob   = c(fit$L[,topics]))
  dat$topic <- factor(dat$topic,topics)
  return(dat)
}
