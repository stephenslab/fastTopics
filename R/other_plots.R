#' @title Plot Progress of Model Fitting Over Time
#'
#' @description Create a plot showing improvement in one or more
#'   Poisson NMF or multinomial topic model fits over time.
#'
#' @details The horizontal axis shows the recorded runtime (in s), and
#' the vertical axis shows some quantity measuring the quality of the
#' fit: the log-likelihood, deviance or maximum residual of the
#' Karush-Kuhn-Tucker (KKT) first-order optimality conditions. To
#' better visualize log-likelihoods and deviances, log-likelihood and
#' deviance differences are shown on the logarithmic scale.
#' Differences are calculated with respect to the best value achieved
#' over all the fits compared.
#'
#' Note that only minimal argument checking is performed.
#' 
#' @param fits An object of class \code{"poisson_nmf_fit"} or
#'   \code{"multinom_topic_model_fit"}, or a non-empty, named list in
#'   which each all list elements are objects of class
#'   \code{"poisson_nmf_fit"} or all objects of class
#'   \code{"multinom_topic_model_fit"}.
#'
#' @param x Choose \code{"timing"} to plot improvement in the solution
#'   over time, or choose \code{"iter"} to plot improvement in the
#'   solution per iteration.
#' 
#' @param y Column of the "progress" data frame used to assess
#'   progress of the Poisson NMF optimization method(s). Should be one
#'   of \code{"loglik"} (Poisson NMF or multinomial topic model
#'   log-likelihood), \code{"dev"} (deviance) or \code{"res"} (maximum
#'   residual of KKT conditions). The deviance is only valid for Poisson
#'   NMF model fits.
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
#'   If fewer line types than \dQuote{fits} are given, the line types are
#'   recycled.
#'
#' @param linesizes Line sizes used to draw progress curves; passed as
#'   the \code{values} input to \code{\link[ggplot2]{scale_size_manual}}.
#'   If fewer line sizes than \dQuote{fits} are given, the line sizes are
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
#' @param theme The \sQuote{ggplot2} \dQuote{theme}.
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
plot_progress <-
  function (fits, x = c("timing","iter"), y = c("loglik","dev","res"),
            add.point.every = 20,
            colors = c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2",
                       "#D55E00","#CC79A7"),
            linetypes = "solid", linesizes = 0.5, shapes = 19, fills = "white",
            e = 0.01, theme = function() theme_cowplot(12)) {

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
                 "or a non-empty, named list in which all list elements are",
                 "objects of class \"poisson_nmf_fit\" or all of class",
                 "\"multinom_topic_model_fit\"")
    if (!(is.list(fits) & !is.null(names(fits)) & length(fits) > 0))
      stop(msg)
    if (!(all(sapply(fits,function(x)inherits(x,"poisson_nmf_fit"))) |
          all(sapply(fits,function(x)inherits(x,"multinom_topic_model_fit")))))
      stop(msg)
    if (!all(nchar(names(fits)) > 0))
      stop(msg)
  }
      
  # Check and process input arguments "x" and "y".
  x <- match.arg(x)
  y <- match.arg(y)
  if (y == "dev" & !inherits(fits[[1]],"poisson_nmf_fit"))
    stop("y = \"dev\" is only valid for Poisson NMF model fits")
  if (y == "loglik" & inherits(fits[[1]],"multinom_topic_model_fit"))
    y <- "loglik.multinom"
 
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

  # Combine the progress information from all the fits into one data
  # frame.
  pdat <- prepare_progress_plot_data(fits,e)

  # Create the plot showing the improvement in the log-likelihood,
  # deviance, or maximum KKT residual over time.
  return(create_progress_plot(pdat,x,y,add.point.every,colors,linetypes,
                              linesizes,shapes,fills,theme))
}

# Used by plot_progress to create a data frame suitable for plotting
# with 'ggplot'.
prepare_progress_plot_data <- function (fits, e) {
  n     <- length(fits)
  labels <- names(fits)
  for (i in 1:n) {
    y         <- fits[[i]]$progress
    y         <- cbind(data.frame("method" = labels[i]),y)
    y$timing  <- cumsum(y$timing)
    fits[[i]] <- y
  }
  out                 <- do.call(rbind,fits)
  out$method          <- factor(out$method,labels)
  out$loglik          <- max(out$loglik) - out$loglik + e
  out$loglik.multinom <- max(out$loglik.multinom) - out$loglik.multinom + e
  out$dev             <- out$dev - min(out$dev) + e
  return(out)
}

# Used by plot_progress to create the plot.
create_progress_plot <- function (pdat, x, y, add.point.every, colors,
                                  linetypes, linesizes, shapes, fills,
                                  theme) {
  rows <- which(pdat$iter %% add.point.every == 1)
  if (x == "timing")
    xlab <- "runtime (s)"
  else if (x == "iter")
    xlab <- "iteration"
  if (y == "res")
    ylab <- "max KKT residual"
  else if (y == "dev")
    ylab <- paste("distance from best deviance")
  else if (y == "loglik" | y == "loglik.multinom")
    ylab <- paste("distance from best loglik")
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
#' @title Plot Log-Likelihood Versus Rank 
#'
#' @description Create a plot showing the improvement in the
#'   log-likelihood as the rank of the matrix factorization or the
#'   number of topics (\dQuote{k}) increases.
#' 
#' @param fits A list with 2 more list elements, in which each list
#'   element is an object of class \code{"poisson_nmf_fit"} or
#'   \code{"multinom_topic_model_fit"}. If two or more fits share the
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
  msg <- paste("Input argument \"fits\" should be a list of length 2 or more",
               "in which all list elements are Poisson NMF fits or all",
               "multinomial topic model fits")
  if (!(is.list(fits) & length(fits) > 1))
    stop(msg)
  if (!(all(sapply(fits,function (x) inherits(x,"poisson_nmf_fit"))) |
        all(sapply(fits,function (x) inherits(x,"multinom_topic_model_fit")))))
    stop(msg)
  n <- length(fits)
  names(fits) <- paste0("fit",1:n)
  dat   <- compare_fits(fits)[c("k","loglik.diff")]
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
#'   settings; most of the \dQuote{heavy lifting} is done by
#'   \sQuote{ggplot2} (specifically, function
#'   \code{\link[ggplot2]{geom_boxplot}} in the \sQuote{ggplot2}
#'   package). For more control over the plot's appearance, the plot can
#'   be customized by modifying the \code{ggplot_call} and
#'   \code{plot_grid_call} arguments.
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

