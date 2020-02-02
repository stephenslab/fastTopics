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
#' @param color Describe "color" input argument here.
#'
#' @param linetype Describe "linetype" input argument here.
#'
#' @param e Describe "e" input argument here.
#'
#' @return Describe return value here.
#'
#' @seealso \code{\link{fit_poisson_nmf}}
#' 
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom cowplot theme_cowplot
#' 
#' @export
#' 
plot_progress_poisson_nmf <- function (fits, color, linetype, size, e = 0.01) {

  # Check input "fits".
  if (!(is.list(fits) & !is.null(names(fits))))
    stop("Input argument \"fits\" should be a list with named elements")
  if (!all(sapply(fits,function (x) inherits(x,"poisson_nmf_fit"))))
    stop("Input argument \"fit\" should be a list in which each list ",
         "element is an object of class \"poisson_nmf_fit\"")
  if (!all(nchar(names(fits)) > 0))
    stop("All names for elements of \"fit\" should be non-empty")
  
  # Prepare the data for plotting.
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
                    dev    = dev - min(dev) + e,
                    delta  = pmax(delta.f,delta.l))

  # Create the plot showing the improvement in the log-likelihood over
  # time.
  p1 <- ggplot(pdat,aes_string(x = "timing",y = "loglik",color = "method")) +
    geom_line(size = 1) +
    scale_y_continuous(trans = "log10") +
    theme_cowplot()

  # Create the plot showing the improvement in the deviance over time.
  # TO DO.

  # Create the plot showing the evolution in the KKT residual over time.
  # TO DO.

  # Create the plot showing the maximum change in the factors and
  # loadings over time.
  return(p1)
}
