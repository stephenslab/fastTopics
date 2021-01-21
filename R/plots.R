#' @title Plot Progress of Poisson NMF Optimization Over Time
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
#' @param fits An object of class \code{"poisson_nmf_fit"} or
#'   \code{"multinom_topic_model_fit"}, or a non-empty, named list in
#'   which each list element is an object of class
#'   \code{"poisson_nmf_fit"} or \code{"multinom_topic_model_fit"}.
#'   Multinomial topic model fits are automatically converted to
#'   equivalent Poisson NMF fits using \code{\link{multinom2poisson}}.
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
#'   \dQuote{fits} are given, the shapes are recycled.
#' 
#' @param fills Fill colours used to draw points at the selected
#'   iterations; passed as the \code{values} input to
#'   \code{\link[ggplot2]{scale_fill_manual}}. If fewer fill colours
#'   than \dQuote{fits} are given, the fill colours are recycled.
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
            e = 0.01, theme = function() theme_cowplot(12)) {

  # CHECK & PROCESS INPUTS
  # ----------------------
  # Check and process input "fits". It should either be an object of
  # class poisson_nmf_fit, or a list of poisson_nmf_fit objects.
 if (inherits(fits,"poisson_nmf_fit") |
     inherits(fits,"multinom_topic_model_fit")) {
    fit.name    <- deparse(substitute(fits))          
    fits        <- list(fits)
    names(fits) <- fit.name
  } else {
    msg <- paste("Input argument \"fits\" should either be an object of",
                 "class \"poisson_nmf_fit\" or \"multinom_topic_model_fit\"",
                 "or a non-empty, named list in which each list element is",
                 "an object of class \"poisson_nmf_fit\" or",
                 "\"multinom_topic_model_fit\"")
    if (!(is.list(fits) & !is.null(names(fits)) & length(fits) > 0))
      stop(msg)
    if (!all(sapply(fits,function (x) inherits(x,"poisson_nmf_fit") |
                                      inherits(x,"multinom_topic_model_fit"))))
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

#' @rdname plot_loglik_vs_rank
#' 
#' @title Plot Log-Likelihood Versus Poisson NMF Rank 
#'
#' @description Create a plot showing the improvement in the Poisson
#'   NMF log-likelihood as the rank of the matrix factorization or the
#'   number of topics (\dQuote{k}) increases.
#' 
#' @param fits A list with 2 more list elements, in which each list
#'   element is an object of class \code{"poisson_nmf_fit"} or
#'   \code{"multinom_topic_model_fit"}. If two or more fits shares the
#'   same rank, or number of topics, the largest log-likelihood is
#'   plotted.
#'
#' @param ggplot_call The function used to create the plot. Replace
#'   \code{loglik_vs_rank_ggplot_call} with your own function to
#'   customize the appearance of the plot.
#'
#' @return A \code{ggplot} object.
#'
#' @export
#'
plot_loglik_vs_rank <- function (fits,
                                 ggplot_call = loglik_vs_rank_ggplot_call) {
  msg <- paste("Input argument \"fits\" should be a list of length 2 or more ",
               "in which each list element is an object of class",
               "\"poisson_nmf_fit\" or \"multinom_topic_model_fit\"")
  if (!(is.list(fits) & length(fits) > 1))
    stop(msg)
  if (!all(sapply(fits,function (x)
                         inherits(x,"poisson_nmf_fit") |
                         inherits(x,"multinom_topic_model_fit"))))
    stop(msg)
  n <- length(fits)
  names(fits) <- paste0("fit",1:n)
  for (i in 1:n)
    if (inherits(fits[[i]],"multinom_topic_model_fit"))
      fits[[i]] <- multinom2poisson(fits[[i]])
  dat   <- compare_poisson_nmf_fits(fits)[c("k","loglik.diff")]
  dat$k <- factor(dat$k)
  y     <- tapply(dat$loglik.diff,dat$k,max)
  dat   <- data.frame(x = as.numeric(names(y)),y = y)
  return(loglik_vs_rank_ggplot_call(dat))
}

#' @rdname plot_loglik_vs_rank
#'
#' @param dat A data frame passed as input to
#'   \code{\link[ggplot2]{ggplot}}, containing, at a minimum, columns
#'   \dQuote{x} and \dQuote{y}.
#'
#' @param font.size Font size used in plot.
#' 
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 scale_x_continuous
#' @importFrom ggplot2 labs
#' @importFrom cowplot theme_cowplot
#' 
#' @export
#' 
loglik_vs_rank_ggplot_call <- function (dat, font.size = 9)
  return(ggplot(dat,aes_string(x = "x",y = "y")) +
         geom_line() +
         geom_point() +
         scale_x_continuous(breaks = dat$x) +
         labs(x = "rank, or number of topics (k)",
              y = "log-likelihood difference") +
         theme_cowplot(font.size))

#' @rdname loadings_plot
#' 
#' @title Loadings Plot
#'
#' @description Generate one or more barcharts to visualize the
#'   relationship between the loadings or mixture proportions and a
#'   selected categorical variable (a factor).
#'
#' @details This is a lightweight interface primarily intended to
#'   expedite creation of boxplots for investigating relationships
#'   between topics and a categorical variables of interest without
#'   having to spend a great deal of time worrying about the plotting
#'   settings; most of the "heavy lifting" is done by ggplot2
#'   (specifically, function \code{\link[ggplot2]{geom_boxplot}} in the
#'   ggplot2 package). For more control over the plot's appearance, the
#'   plot can be customized by modifying the \code{ggplot_call}
#'   and \code{plot_grid_call} arguments.
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
         "\"multinom_topic_model_fit\"")
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
      plots[[i]] <- loadings_plot(fit,x,k[i],ggplot_call,NULL)
    return(plot_grid_call(plots))
  }
}

#' @rdname loadings_plot
#'
#' @param dat A data frame passed as input to
#'   \code{\link[ggplot2]{ggplot}}, containing, at a minimum, columns
#'   \dQuote{x} and \dQuote{loading}.
#'
#' @param topic.label The name or number of the topic being plotted.
#'   Only used to determine the plot title.
#' 
#' @param font.size Font size used in plot.
#' 
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 geom_boxplot
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_text
#' @importFrom cowplot theme_cowplot
#' 
#' @export
#' 
loadings_plot_ggplot_call <- function (dat, topic.label, font.size = 9)
  ggplot(dat,aes_string(x = "x",y = "loading")) +
    geom_boxplot(width = 0.25,size = 0.4,outlier.shape = NA) +
    labs(x = "",y = "loading",title = paste("topic",topic.label)) +
    theme_cowplot(font.size) +
    theme(axis.text.x = element_text(angle = 45,hjust = 1),
          plot.title  = element_text(size = font.size,face = "plain"))

#' @rdname volcano_plot
#'
#' @title Volcano Plot
#'
#' @description Create one or more "volcano" plots to visualize the
#'   results of a differential count analysis using a topic model. A
#'   volcano plot is a scatterplot in which the log-fold change,
#'   estimated using a multinomial topic model, is plotted against the
#'   p-value or z-score. 
#'
#' @details The colour of the points is varied by the average count,
#' on the logarithmic scale; since the evidence (z-score or p-value)
#' typically increases with more observed counts, the variables with
#' smallest average counts should usually appear toward the bottom of
#' the volcano plot. Only points above a specified z-score (or
#' p-value) quantile are labeled. Note that points with an average
#' count of zero or less are not shown.
#'
#' To better accommodate situations in which some z-scores (or
#' p-values) are much larger than all the others, the z-scores and
#' -log10 p-values are plotted on the square-root scale. To change
#' this, as well as other aspects, of the plot, replace
#' \code{volcano_plot_ggplot_call} with your own function; see input
#' argument \dQuote{ggplot_call}.
#'
#' The \dQuote{ggrepel} package is used to arrange the labels in a
#' visually attractive manner.
#'
#' Use interactive volcano plot is created using the \dQuote{plotly}
#' package. The "hover text" shows the label (see input argument
#' \dQuote{labels}) and detailed log-fold change statistics as they
#' were calculated by \code{\link{diff_count_analysis}}.
#' 
#' @param diff_count_result An object of class
#'   \dQuote{topic_model_diff_count}, usually an output from
#'   \code{\link{diff_count_analysis}}.
#'
#' @param k The topic, or topics, selected by number or name. When not
#'   specified, all topics are plotted.
#'
#' @param labels Character vector specifying how the points in the
#'   volcano plot are labeled. This should be a character vector with
#'   one entry per log-fold change statistic (row of
#'   \code{diff_count_result$beta}). When not specified, the row names
#'   of \code{diff_count_result$beta} are used, if available. Labels are
#'   added to the plot using \code{\link[ggrepel]{geom_text_repel}}.
#'
#' @param y The measure of support to plot in the y-axis: use \code{y
#'   = "zscore"} to plot the z-score magnitudes; use \code{y = "pvalue"}
#'   to plot the -log1 p-values.
#' 
#' @param betamax Truncate the log-fold change statistics
#'   (\code{beta}) by this amount. Any statistics greater in magnitude
#'   than \code{betamax} are set to \code{betamax} or \code{-betamax}.
#'
#' @param label_above_quantile Only z-scores or p-values (depending on
#'   choice of \code{y}) above this quantile are labeled in the volcano
#'   plot. \code{\link[ggrepel]{geom_text_repel}} will attempt to label
#'   all points when \code{label_above_quantile = 0}. When
#'   \code{label_above_quantile = Inf}, no points are labeled.
#'
#' @param subsample_below_quantile A number between 0 and 1. If
#'   greater than zero, log-fold change statistics with z-scores or
#'   p-values below this quantile will be subsampled according to
#'   \code{subsample_rate}. This is useful for large data sets to to
#'   reduce the number of points plotted.
#'
#' @param subsample_rate A number between 0 and 1 giving the
#'   proportion of log-fold change statistics with "small" z-scores that
#'   are included in the plot, uniformly at random. This is only used if
#'   \code{subsample_below_quantile} is greater than zero.
#'
#' @param max.overlaps Argument passed to
#'   \code{\link[ggrepel]{geom_text_repel}}.
#' 
#' @param ggplot_call The function used to create the plot. Replace
#'   \code{volcano_plot_ggplot_call} with your own function to customize
#'   the appearance of the plot.
#'
#' @param plot_grid_call When multiple topics are selected, this is
#'   the function used to arrange the plots into a grid using
#'   \code{\link[cowplot]{plot_grid}}. It should be a function accepting
#'   a single argument, \code{plots}, a list of \code{ggplot} objects.
#'
#' @return A \code{ggplot} object or a \code{plotly} object.
#'
#' @importFrom cowplot plot_grid
#'
#' @export
#' 
volcano_plot <-
  function (diff_count_result, k, labels, y = c("zscore", "pvalue"),
            betamax = 10, label_above_quantile = 0.99,
            subsample_below_quantile = 0, subsample_rate = 0.1,
            max.overlaps = Inf, ggplot_call = volcano_plot_ggplot_call,
            plot_grid_call = function (plots) do.call(plot_grid,plots)) {
    
  # Check and process input arguments.
  y <- match.arg(y)
  if (!inherits(diff_count_result,"topic_model_diff_count"))
    stop("Input \"diff_count_result\" should be an object of class ",
         "\"topic_model_diff_count\"")
  beta <- diff_count_result$beta
  if (missing(k))
    k <- seq(1,ncol(beta))
  if (missing(labels)) {
    if (!is.null(rownames(beta)))
      labels <- rownames(beta)
    else
      labels <- as.character(seq(1,nrow(beta)))
  }
  if (!(is.character(labels) & length(labels) == nrow(beta)))
    stop("Input argument \"labels\", when specified, should be a character ",
         "vector with one entry per log-fold change statistic (column of ",
         "the counts matrix)")

  if (length(k) == 1) {
    if (y == "zscore")
      y.label <- "|z-score|"
    else if (y == "pvalue")
      y.label <- "-log10 p-value"
    dat <- compile_volcano_plot_data(diff_count_result,k,labels,y,betamax,
                                     label_above_quantile,
                                     subsample_below_quantile,subsample_rate)
    return(ggplot_call(dat,y.label,k,max.overlaps))
  } else {

    # Create a volcano plot for each selected topic, and combine them
    # using plot_grid. This is done by recursively calling volcano_plot.
    m     <- length(k)
    plots <- vector("list",m)
    names(plots) <- k
    for (i in 1:m)
      plots[[i]] <- volcano_plot(diff_count_result,k[i],labels,y,betamax,
                                 label_above_quantile,subsample_below_quantile,
                                 subsample_rate,max.overlaps,ggplot_call,NULL)
    return(plot_grid_call(plots))
  }
}

#' @rdname volcano_plot
#'
#' @param x An object of class \dQuote{topic_model_diff_count},
#'   usually an output from \code{\link{diff_count_analysis}}.
#'
#' @param \dots Additional arguments passed to \code{volcano_plot}.
#' 
#' @importFrom graphics plot
#' 
#' @method plot topic_model_diff_count
#'
#' @export
#'
plot.topic_model_diff_count <- function (x, ...)
  volcano_plot(x,...)

#' @rdname volcano_plot
#'
#' @param file Save the interactive volcano plot to this HTML
#'   file using \code{\link[htmlwidgets]{saveWidget}}.
#' 
#' @param width Width of the plot in pixels. Passed as argument
#'   \dQuote{width} to \code{\link[plotly]{plot_ly}}.
#'
#' @param height Height of the plot in pixels. Passed as argument
#'   \dQuote{height} to \code{\link[plotly]{plot_ly}}.
#'
#' @param title The text used for the plot title.
#'
#' @param plot_ly_call The function used to create the plot. Replace
#'   \code{volcano_plot_ly_call} with your own function to customize
#'   the appearance of the interactive plot.
#'
#' @importFrom htmlwidgets saveWidget
#' 
#' @export
#' 
volcano_plotly <- function (diff_count_result, k, file, labels,
                            y = c("zscore", "pvalue"), betamax = 10,
                            subsample_below_quantile = 0,
                            subsample_rate = 0.1, width = 600, height = 500,
                            title = paste("topic",k),
                            plot_ly_call = volcano_plot_ly_call) {

  # Check and process input arguments.
  y <- match.arg(y)
  if (!inherits(diff_count_result,"topic_model_diff_count"))
    stop("Input \"diff_count_result\" should be an object of class ",
         "\"topic_model_diff_count\"")
  beta <- diff_count_result$beta
  if (missing(labels)) {
    if (!is.null(rownames(beta)))
      labels <- rownames(beta)
    else
      labels <- as.character(seq(1,nrow(beta)))
  }
  if (!(is.character(labels) & length(labels) == nrow(beta)))
    stop("Input argument \"labels\", when specified, should be a character ",
         "vector with one entry per log-fold change statistic (column of ",
         "the counts matrix)")

  # Compile the plotting data.
  dat <- compile_volcano_plot_data(diff_count_result,k,labels,y,betamax,
                                   0,subsample_below_quantile,subsample_rate)

  # Create the interactive volcano plot using plotly.
  if (y == "zscore")
    y.label <- "|z-score|"
  else if (y == "pvalue")
    y.label <- "-log10 p-value"
  p <- volcano_plot_ly_call(dat,y.label,title,width,height)
  if (!missing(file))
    saveWidget(p,file,selfcontained = TRUE,title = title)
  return(p)
}

#' @rdname volcano_plot
#'
#' @param dat A data frame passed as input to
#'   \code{\link[ggplot2]{ggplot}}, containing, at a minimum, columns
#'   \dQuote{label}, \dQuote{mean}, \dQuote{beta}, \dQuote{se},
#'   \dQuote{z}, \dQuote{pval} and \dQuote{label}.
#' 
#' @param y.label Label to use in the plot for \dQuote{dat$y}.
#' 
#' @param topic.label The name or number of the topic being plotted.
#'   Only used to determine the plot title.
#' 
#' @param max.overlaps Argument passed to
#'   \code{\link[ggrepel]{geom_text_repel}}.
#' 
#' @param font.size Font size used in plot.
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 scale_x_continuous
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom ggplot2 scale_fill_gradient2
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 expansion
#' @importFrom ggplot2 element_text
#' @importFrom ggrepel geom_text_repel
#' @importFrom cowplot theme_cowplot
#' 
#' @export
#' 
volcano_plot_ggplot_call <- function (dat, y.label, topic.label,
                                      max.overlaps = Inf, font.size = 9)
  ggplot(dat,aes_string(x = "beta",y = "y",fill = "mean",label = "label")) +
    geom_point(color = "white",stroke = 0.3,shape = 21,na.rm = TRUE) +
    scale_x_continuous(expand = expansion(mult = 0.2)) +
    scale_y_continuous(trans = "sqrt",
      breaks = c(0,1,2,5,10,20,50,100,200,500,1e3,2e3,5e3,1e4,2e4,5e4)) +
    scale_fill_gradient2(low = "deepskyblue",mid = "gold",high = "orangered",
                         na.value = "gainsboro",
                         midpoint = mean(range(dat$mean))) +
    geom_text_repel(color = "black",size = 2.25,fontface = "italic",
                    segment.color = "black",segment.size = 0.25,
                    min.segment.length = 0,max.overlaps = max.overlaps,
                    na.rm = TRUE) +
    labs(x = "log-fold change (\u03b2)",y = y.label,fill = "log10 mean",
         title = paste("topic",topic.label)) +
    theme_cowplot(font.size) +
    theme(plot.title = element_text(size = font.size,face = "plain"))   

#' @rdname volcano_plot
#'
#' @importFrom plotly plot_ly
#' @importFrom plotly hide_colorbar
#' @importFrom plotly layout
#' 
#' @export
#' 
volcano_plot_ly_call <- function (dat, y.label, title, width, height) {
  p <- plot_ly(data = dat,x = ~beta,y = ~sqrt(y),color = ~mean,
               colors = c("deepskyblue","gold","orangered"),
               text = ~sprintf(paste0("%s\nmean: %0.3f\n\u03b2: %+0.3f\n",
                                      "s.e.: %0.3f\nz: %+0.3f\n-log10p: %0.2f"),
                               label,10^mean,beta,se,z,pval),
               type = "scatter",mode = "markers",hoverinfo = "text",
               width = width,height = height,
               marker = list(line = list(color = "white",width = 1),size = 7.5))
  p <- hide_colorbar(p)
  p <- layout(p,xaxis = list(title = "log-fold change (\u03b2)",
                             zeroline = FALSE,showgrid = FALSE),
              yaxis = list(title = paste("sqrt",y.label),
                           zeroline = FALSE,showgrid = FALSE),
              hoverlabel = list(bgcolor = "white",bordercolor = "black",
                                font = list(color = "black",family = "arial",
                                            size = 12)),
              font = list(family = "arial",size = 12),
              showlegend = FALSE,title = title)
  return(p)
}

# This is used by volcano_plot and volcano_plotly to compile the data
# frame passed to ggplot.
#
#' @importFrom stats quantile
compile_volcano_plot_data <- function (diff_count_result, k, labels, y,
                                       betamax, label_above_quantile,
                                       subsample_below_quantile,
                                       subsample_rate) {
  dat <- with(diff_count_result,
                data.frame(label = labels,
                           mean  = colmeans,
                           beta  = beta[,k],
                           se    = se[,k],
                           z     = Z[,k],
                           pval  = pval[,k],
                           y     = 0,
                           stringsAsFactors = FALSE))
  if (y == "zscore")
    dat$y <- abs(diff_count_result$Z[,k])
  else if (y == "pvalue")
    dat$y <- diff_count_result$pval[,k]
  rows     <- which(dat$mean > 0)
  dat      <- dat[rows,]
  dat$mean <- log10(dat$mean)
  dat <- transform(dat,beta = sign(beta) * pmin(betamax,abs(beta)))
  if (is.infinite(label_above_quantile))
    y0 <- Inf
  else
    y0 <- quantile(dat$y,label_above_quantile)
  dat$label[dat$y < y0] <- ""
  if (subsample_below_quantile > 0) {
    y0    <- quantile(dat$y,subsample_below_quantile)
    rows1 <- which(dat$y >= y0)
    rows2 <- which(dat$y < y0)
    rows2 <- sample(rows2,ceiling(subsample_rate * length(rows2)))
    rows  <- sort(c(rows1,rows2))
    message(sprintf("%d out of %d data points will be included in plot",
                    length(rows),nrow(dat)))
    dat   <- dat[rows,]
  }
  return(dat)
}

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
#' @importFrom stats prcomp
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

#' @title t-SNE from Poisson NMF or Multinomial Topic Model
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
#'   importance, or "weight", of the respective topic in computing the
#'   embedding. When \code{scaling} is \code{NULL}, no re-scaling is
#'   performed. Note that this scaling will have no effect if
#'   \code{normalize = TRUE}.
#' 
#' @param pca Whether to perform a PCA processing stepe in t-SNE;
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
    
  # Check and process input arguments.
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

#' @title t-SNE Plot
#'
#' @description Visualize the "structure" of the Poisson NMF loadings
#'   or multinomial topic model mixture proportions by projection onto a
#'   2-d surface. Samples in the projection are colored according to the
#'   their loadings/proportions. By default, t-SNE is used to compute
#'   the 2-d embedding from the loadings or mixture proportions.
#'
#' @details This is a lightweight interface primarily intended to
#'   expedite creation of scatterplots for visualizing the loadings or
#'   mixture proportions in 2-d; most of the "heavy lifting" is done by
#'   ggplot2. The 2-d embedding itself is computed by invoking function
#'   \code{\link{tsne_from_topics}} (unless the "tsne" input is
#'   provided). For more control over the plot's appearance, the plot
#'   can be customized by modifying the \code{ggplot_call} and
#'   \code{plot_grid_call} arguments.
#'
#'   An effective 2-d visualization may also necessitate some fine-tunning
#'   of the t-SNE settings, such as the "perplexity", or the number of
#'   samples included in the plot. The t-SNE settings can be controlled
#'   by the additional arguments (\dots) passed to
#'   \code{tsne_from_topics}; see \code{\link{tsne_from_topics}} for
#'   details. Alternatively, a 2-d embedding may be pre-computed, and
#'   passed as argument \code{tsne} to \code{tsne_plot}.
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
#'   "d1", "d2" (the first and second dimensions in the 2-d embedding),
#'   and "loading".
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

#' @title Structure Plot
#'
#' @description Create a "Structure plot" from a multinomial topic
#'   model fit. The Structure plot represents the estimated mixture
#'   proportions of each sample as a stacked barchart, with bars of
#'   different colors representing different topics. Consequently,
#'   samples that have similar mixture proportions have similar amounts
#'   of each color.
#'
#' @details The name "Structure plot" comes from its widespread use in
#'   population genetics to visualize the results of the Structure
#'   software (Rosenberg \emph{et al}, 2002).
#' 
#'   For most uses of the Structure plot in population genetics, there
#'   is usually some grouping of the samples (e.g., assignment to
#'   pre-defined populations) that guides arrangement of the samples
#'   along the horizontal axis in the bar chart. In other applications,
#'   such as analysis of gene expression data, no pre-defined grouping
#'   exists. Therefore, a "smart" arrangement of the samples is, by
#'   default, generated automatically by performing a 1-d t-SNE
#'   embedding of the samples.
#'
#'   Alternatively, a categorical variable---the grouping---may be
#'   provided, in which case the samples are arranged according to that
#'   grouping, then arranged within each group using t-SNE.
#'
#' @param fit An object of class \dQuote{poisson_nmf_fit} or
#'   \dQuote{multinom_topic_model_fit}. If a Poisson NMF fit is provided
#'   as input, the corresponding multinomial topic model fit is
#'   automatically recovered using \code{\link{poisson2multinom}}.
#'
#' @param n The maximum number of samples (rows of the loadings
#'   matrix) to include in the plot. Including large numbers (e.g.,
#'   thousands) of samples is not recommended because it dramatically
#'   slows down the t-SNE computation, and typically there is little to
#'   no benefit in including large number of samples in the plot due to
#'   screen resolution limits. Ignored if \code{rows} is provided.
#'
#' @param rows Ordering of rows (samples) in plot, after they have
#'   been grouped. Generated automatically from 1-d t-SNE if not
#'   provided.
#'
#' @param grouping Optional categorical variable (factor) with one
#'   entry for each row of the loadings matrix \code{fit$L} defining a
#'   grouping of the samples (rows). The samples (rows) are arranged
#'   along the horizontal axis according to this grouping, then within
#'   each group according to \code{rows} or, if not provided, according
#'   to the 1-d t-SNE embedding.
#' 
#' @param topics Top-to-bottom ordering of the topics in the structure
#'   plot; topics[1] is shown on the top, topics[2] is shown next, and
#'   so on. If the ordering of the topics is not specified, the topics
#'   are automaticcally ordered so that the topics with the greatest
#'   "mass" are at shown at the bottom of the plot.
#' 
#' @param colors Colors used to draw topics in Structure plot:
#'   \code{colors[1]} is the colour used to draw \code{topics[1]},
#'   \code{colors[2]} is the colour used to draw \code{topics[2]}, and
#'   so on. The default colour setting is the from
#'   \url{https://colorbrewer2.org} (qualitative data, "9-class Set1").
#'
#' @param gap The horizontal spacing between groups. Ignored if
#'   \code{grouping} is not provided.
#'
#' @param perplexity t-SNE perplexity parameter, passed as argument
#'   \dQuote{perplexity} to \code{\link{tsne_from_topics}}. When
#'   code{grouping} is provided, \code{perplexity} may be a vector of
#'   settings of length equal to \code{levels(grouping)}, which allows
#'   for a different t-SNE perplexity setting for each group.
#'
#' @param eta t-SNE learning rate parameter, passed as argument
#'   \dQuote{eta} to \code{\link{tsne_from_topics}}. When
#'   code{grouping} is provided, \code{eta} may be a vector of
#'   settings of length equal to \code{levels(grouping)}, which allows
#'   for a different t-SNE learning rate for each group.
#' 
#' @param ggplot_call The function used to create the plot. Replace
#'   \code{structure_plot_ggplot_call} with your own function to
#'   customize the appearance of the plot.
#'
#' @param \dots Additional arguments passed to \code{structure_plot}
#'   (for the \code{plot} method) or \code{\link{tsne_from_topics}} (for
#'   function \code{structure_plot}).
#'
#' @return A \code{ggplot} object.
#'
#' @references
#'
#' Dey, K. K., Hsiao, C. J. and Stephens, M. (2017). Visualizing the
#' structure of RNA-seq expression data using grade of membership
#' models. \emph{PLoS Genetics} \bold{13}, e1006599.
#'
#' Rosenberg, N. A., Pritchard, J. K., Weber, J. L., Cann, H. M.,
#' Kidd, K. K., Zhivotovsky, L. A. and Feldman, M. W. (2002). Genetic
#' structure of human populations. \emph{Science} \bold{298},
#' 2381–2385.
#'
#' @export
#'
structure_plot <-
  function (fit, n = 2000, rows, grouping, topics,
            colors = c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00",
                       "#ffff33","#a65628","#f781bf","#999999"),
            gap = 1, perplexity = 100, eta = 200,
            ggplot_call = structure_plot_ggplot_call, ...) {

  # Check and process inputs.
  if (!(inherits(fit,"poisson_nmf_fit") |
        inherits(fit,"multinom_topic_model_fit")))
    stop("Input \"fit\" should be an object of class \"poisson_nmf_fit\" or ",
         "\"multinom_topic_model_fit\"")
  if (inherits(fit,"poisson_nmf_fit"))
    fit <- poisson2multinom(fit)
  n0 <- nrow(fit$L)
  k  <- ncol(fit$L)
  if (missing(grouping))
    grouping <- factor(rep(1,n0))
  groups     <- levels(grouping)
  num_groups <- length(groups)
  if (missing(topics))
    topics <- order(colMeans(fit$L))
  if (length(perplexity) == 1)
    perplexity <- rep(perplexity,num_groups)
  if (length(eta) == 1)
    eta <- rep(eta,num_groups)
  names(perplexity) <- groups
  names(eta)        <- groups
  if (num_groups == 1) {

    # If the ordering of the rows is not provided, determine an
    # ordering by computing a 1-d embedding from the mixture
    # proportions.
    if (missing(rows)) {
      out  <- tsne_from_topics(fit,1,n,perplexity = perplexity,eta = eta,...)
      rows <- out$rows[order(out$Y)]
    }

    # Prepare the data for plotting. Note that it is important to
    # subset the rows *after* recovering the parameters of the
    # multinomial topic model.
    dat <- compile_structure_plot_data(fit$L[rows,],topics)

    # Create the Structure plot.
    return(ggplot_call(dat,colors))
  } else {

    # If the ordering of the rows is not provided, determine an
    # ordering by computing a 1-d embedding from the mixture
    # proportions, separately for each group.
    if (missing(rows)) {
      rows <- NULL
      n    <- round(n/n0*table(grouping))
      for (j in levels(grouping)) {
        i    <- which(grouping == j)
        out  <- tsne_from_topics(select_loadings(fit,i),1,n[j],
                                 perplexity = perplexity[j],
                                 eta = eta[j],...)
        rows <- c(rows,i[out$rows[order(out$Y)]])
      }
    }

    # Prepare the data for plotting. Note that it is important to
    # subset the rows *after* recovering the parameters of the
    # multinomial topic model.
    out <- compile_grouped_structure_plot_data(fit$L[rows,],topics,
                                               grouping[rows],gap)

    # Create the Structure plot.
    return(ggplot_call(out$dat,colors,out$ticks))
  }
}

#' @rdname structure_plot
#'
#' @param x An object of class \dQuote{poisson_nmf_fit} or
#'   \dQuote{multinom_topic_model_fit}. If a Poisson NMF fit is provided
#'   as input, the corresponding multinomial topic model fit is
#'   automatically recovered using \code{\link{poisson2multinom}}.
#'
#' @importFrom graphics plot
#' 
#' @method plot poisson_nmf_fit
#'
#' @export
#'
plot.poisson_nmf_fit <- function (x, ...)
  structure_plot(x,...)

#' @rdname structure_plot
#'
#' @importFrom graphics plot
#' 
#' @method plot multinom_topic_model_fit
#' 
#' @export
#' 
plot.multinom_topic_model_fit <- function (x, ...)
  structure_plot(x,...)

#' @rdname structure_plot
#'
#' @param dat A data frame passed as input to
#'   \code{\link[ggplot2]{ggplot}}, containing, at a minimum, columns
#'   "sample", "topic" and "mixprop": column "sample" contains the
#'   positions of the samples (rows of the loadings matrix) along the
#'   horizontal axis; column "topic" is a topic (corresponding to
#'   columns of the loadings matrix); and column "mixprop" is the mixture
#'   proportion for the given sample.
#'
#' @param ticks The placement of the group labels along the horizontal
#'   axis, and their names. For data that is not grouped, use
#'   \code{ticks = NULL.}
#'
#' @param font.size Font size used in plot.
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
#' @importFrom ggplot2 element_text
#' @importFrom cowplot theme_cowplot
#' 
#' @export
#'
structure_plot_ggplot_call <- function (dat, colors, ticks = NULL,
                                        font.size = 9)
  ggplot(dat,aes_string(x = "sample",y = "mixprop",color = "topic",
                        fill = "topic")) +
    geom_col() +
    scale_x_continuous(limits = c(0,max(dat$sample) + 1),
                       breaks = ticks,labels = names(ticks)) +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors) +
    labs(x = "",y = "mixture proportion") +
    theme_cowplot(font.size) +
    theme(axis.line   = element_blank(),
          axis.ticks  = element_blank(),
          axis.text.x = element_text(angle = 45,hjust = 1))

# This is used by structure_plot to create a data frame suitable for
# plotting with ggplot. Input argument L is the "loadings" matrix from
# a multinomial topic model fit---each row of L is a vector of mixture
# proportions. Input argument "topics" is the vector of the selected
# topics (that is, selected columns of L). The output is a data frame
# with three columns: "sample", a row of L (numeric); "topic", a topic
# (factor); and "mixprop", the mixture proportion for the given sample
# (numeric).
compile_structure_plot_data <- function (L, topics) {
  n   <- nrow(L)
  k   <- length(topics)
  dat <- data.frame(sample  = rep(1:n,times = k),
                    topic   = rep(topics,each = n),
                    mixprop = c(L[,topics]))
  dat$topic <- factor(dat$topic,topics)
  return(dat)
}

# This is used by structure_plot to create a data frame suitable for
# plotting with ggplot when the data are grouped. Input argument L is
# the "loadings" matrix from a multinomiial topic model fit---each row
# of L is a vector of mixture proportions. Input argument "topics" is
# the vector of the selected topics (that is, selected columns of L).
# Input argument "grouping" is a factor with one entry for each row of
# L giving the grouping for the rows of L. The rows of L (and,
# correspondingly, the grouping vector) should be arranged in order of
# the groups; that is, "sort(grouping)" should be the same as
# "grouping". Finally, the value "gap" (zero by default) is added to
# the sample indices in each group to provide a visual spacing of the
# groups.
compile_grouped_structure_plot_data <- function (L, topics, grouping,
                                                 gap = 0) {
  ng    <- nlevels(grouping)
  ticks <- rep(0,ng)
  names(ticks) <- levels(grouping)
  dat   <- NULL
  m     <- 0
  for (j in levels(grouping)) {
    i          <- which(grouping == j)
    out        <- compile_structure_plot_data(L[i,],topics)
    out$sample <- out$sample + m
    n          <- length(i)
    dat        <- rbind(dat,out)
    ticks[j]   <- m + n/2
    m          <- m + n + gap
  }
  return(list(dat = dat,ticks = ticks))
}
