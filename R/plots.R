# NOTES:
#
#  - Caution that very little argument checking is done.
#

#' @title Plot Progress of NMF Optimization Method Over Time
#'
#' @description Add description here.
#'
#' @param fits An object of class \code{"poisson_nmf_fit"}, or a
#'   non-empty, named list in which each list element is an object of
#'   class \code{"poisson_nmf_fit"}.
#'
#' @param y Describe "y" input argument here.
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
#' @param e Describe "e" input argument here.
#'
#' @param theme Describe "theme" input argument here.
#'
#' @return Describe return value here.
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
            linetypes = "solid", linesizes = 0.5, shapes = 21, fills = "white",
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
  labels <- names(fits)
  pdat    <- NULL
  for (i in 1:n) {
    d        <- fits[[i]]$progress
    d        <- cbind(data.frame("method" = labels[i],"iter" = 1:nrow(d)),d)
    d$timing <- cumsum(d$timing)
    pdat     <- rbind(pdat,d)
  }
  pdat$method <- factor(pdat$method,labels)
  pdat$loglik <- max(pdat$loglik) - pdat$loglik + e
  pdat$dev    <- pdat$dev - min(pdat$dev) + e

  # CREATE PLOT
  # -----------
  # Create the plot showing the improvement in the log-likelihood (or
  # deviance) over time.
  y <- "loglik"
  rows <- which(pdat$iter %% add.point.every == 1)
  p1   <- ggplot(pdat,aes_string(x = "timing",y = y,color = "method",
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
    labs(x = "runtime (s)",y = paste("distance from best",y)) +
    theme()
  return(p1)
}
