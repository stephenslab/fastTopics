# NOTES:
#
#  - Caution that very little argument checking is done.
#

#' @title Add title here.
#'
#' @description Add description here.
#'
#' @param fits Describe "fits" input argument here.
#'
#' @param plot.dev Describe "plot.dev" input argument here.
#' 
#' @param color Describe "color" input argument here.
#'
#' @param linetype Describe "linetype" input argument here.
#'
#' @param linesize Describe "linesize" input argument here.
#'
#' @param shape Describe "shape" input argument here.
#' 
#' @param ptsize Describe "ptsize" input argument here.
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
#' @importFrom ggplot2 labs
#' @importFrom cowplot plot_grid
#' @importFrom cowplot theme_cowplot
#' 
#' @export
#' 
plot_progress_poisson_nmf <-
  function (fits,
            plot.dev = FALSE,
            color    = rep(c("#E69F00","#56B4E9","#009E73","#F0E442",
                             "#0072B2","#D55E00","#CC79A7"),
                           length.out = length(fits)),
            linetype = rep("solid",length(fits)),
            linesize = rep(0.75,length(fits)),
            shape    = rep(20,length(fits)),
            ptsize   = 2,
            e        = 0.01,
            theme    = function() theme_cowplot(font_size = 12)) {

  # Check input "fits".
  if (!(is.list(fits) & !is.null(names(fits))))
    stop("Input argument \"fits\" should be a list with named elements")
  if (!all(sapply(fits,function (x) inherits(x,"poisson_nmf_fit"))))
    stop("Input argument \"fit\" should be a list in which each list ",
         "element is an object of class \"poisson_nmf_fit\"")
  if (!all(nchar(names(fits)) > 0))
    stop("All names for elements of \"fit\" should be non-empty")
  
  # Prepare the data for the plots.
  n      <- length(fits)
  labels <- names(fits)
  pdat    <- NULL
  for (i in 1:n) {
    d    <- fits[[i]]$progress
    d    <- cbind(data.frame(method = labels[i],iter = 1:nrow(d)),d)
    d    <- transform(d,timing = cumsum(timing))    
    pdat <- rbind(pdat,d)
  }
  pdat <- transform(pdat,
                    method = factor(method,labels),
                    loglik = max(loglik) - loglik + e,
                    dev    = dev - min(dev) + e)

  # Create the plot showing the improvement in the log-likelihood (or
  # deviance) over time.
  if (plot.dev)
    y <- "dev"
  else
    y <- "loglik"
  p1 <- ggplot(pdat,aes_string(x = "timing",y = y,color = "method",
                               linetype = "method",size = "method")) +
    geom_line() +
    geom_point(data = subset(pdat,iter %% 10 == 1),
               mapping = aes_string(x = "timing",y = y,color = "method",
                                    shape = "method"),
               size = ptsize,inherit.aes = FALSE) +
    scale_y_continuous(trans = "log10") +
    scale_color_manual(values = color) +
    scale_linetype_manual(values = linetype) +
    scale_size_manual(values = linesize) +
    scale_shape_manual(values = shape) +
    labs(x = "runtime (s)",
         y = paste("distance from best",y)) +
    theme()

  # Create the plot showing the evolution in the KKT residual over time.
  p2 <- ggplot(pdat,aes_string(x = "timing",y = "res",color = "method",
                               linetype = "method",size = "method")) +
    geom_line() +
    geom_point(data = subset(pdat,iter %% 10 == 1),
               mapping = aes_string(x = "timing",y = "res",color = "method",
                                    shape = "method"),
               size = ptsize,inherit.aes = FALSE) +
    scale_y_continuous(trans = "log10") +
    scale_color_manual(values = color) +
    scale_linetype_manual(values = linetype) +
    scale_size_manual(values = linesize) +
    scale_shape_manual(values = shape) +
    labs(x = "runtime (s)",
         y = "KKT residual") +
    theme()

  return(plot_grid(p1,p2,labels = c("A","B")))
}
