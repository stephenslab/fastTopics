# NOTES:
#
#  - Caution that very little argument checking is done.
#

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
#' @param fits An object of class \code{"poisson_nmf_fit"}, or a
#'   non-empty, named list in which each list element is an object of
#'   class \code{"poisson_nmf_fit"}.
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
#' @return The return value is a \code{ggplot} object.
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
  function (fits, y = c("loglik","dev","res"), add.point.every = 20,
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

  # Check and process input argument "y".
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
  return(create_progress_plot(pdat,y,add.point.every,colors,linetypes,
                              linesizes,shapes,fills,theme))
}

# Used by plot_progress_poisson_nmf to create a data frame suitable
# for plotting with ggplot.
prepare_progress_plot_data <- function (fits, e) {
  n     <- length(fits)
  labels <- names(fits)
  for (i in 1:n) {
    y         <- fits[[i]]$progress
    y         <- cbind(data.frame("method" = labels[i],"iter" = 1:nrow(y)),y)
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
create_progress_plot <- function (pdat, y, add.point.every, colors, linetypes,
                                  linesizes, shapes, fills, theme) {
  rows <- which(pdat$iter %% add.point.every == 1)
  if (y == "res")
    ylab <- "max KKT residual"
  else
    ylab <- paste("distance from best",y)
  return(ggplot(pdat,aes_string(x = "timing",y = y,color = "method",
                                linetype = "method",size = "method")) +
    geom_line(na.rm = TRUE) +
    geom_point(data = pdat[rows,],
               mapping = aes_string(x = "timing",y = y,color = "method",
                                    fill = "method",shape = "method"),
               inherit.aes = FALSE,na.rm = TRUE) +
    scale_y_continuous(trans = "log10") +
    scale_color_manual(values = colors) +
    scale_linetype_manual(values = linetypes) +
    scale_size_manual(values = linesizes) +
    scale_shape_manual(values = shapes) +
    scale_fill_manual(values = fills) +
    labs(x = "runtime (s)",y = ylab) +
    theme())
}
