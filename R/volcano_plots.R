#' @rdname volcano_plot
#'
#' @title Volcano Plot
#'
#' @description Create one or more "volcano" plots to visualize the
#'   results of a differential count analysis using a topic model. A
#'   volcano plot is a scatterplot in which the log-fold change (LFC),
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
#' package. The \dQuote{hover text} shows the label (see input argument
#' \dQuote{labels}) and detailed LFC statistics as they were
#' calculated by \code{\link{diff_count_analysis}}.
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
#'   one entry per LFC estimate (row of
#'   \code{diff_count_result$beta}). When not specified, the row names
#'   of \code{diff_count_result$beta} are used, if available. Labels are
#'   added to the plot using \code{\link[ggrepel]{geom_text_repel}}.
#'
#' @param y The measure of support to plot in the y-axis: use \code{y
#'   = "zscore"} to plot the z-score magnitudes; use \code{y = "pvalue"}
#'   to plot the -log1 p-values.
#'
#' @param betamax Truncate the LFC statistics (\code{beta}) by this
#'   amount. Any statistics greater than than \code{betamax} are set to
#'   \code{betamax}, and any statistics less than \code{-betamax} are
#'   set to \code{-betamax}.
#'
#' @param label_above_lfc Only z-scores or p-values (depending on
#'   choice of \code{y}) with LFC estimates above this value are
#'   labeled in the volcano plot.
#' 
#' @param label_above_quantile Only z-scores or p-values (depending on
#'   choice of \code{y}) above this quantile are labeled in the volcano
#'   plot.
#'
#' @param subsample_below_quantile A number between 0 and 1. If
#'   greater than zero, LFC estimates with z-scores or p-values below
#'   this quantile will be subsampled according to
#'   \code{subsample_rate}. This can be helpful to reduce the number of
#'   points plotted for a large data set.
#'
#' @param subsample_rate A number between 0 and 1 giving the
#'   proportion of LFC estimates with small z-scores that are included
#'   in the plot, uniformly at random. This is only used if
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
            betamax = 15, label_above_lfc = 0, label_above_quantile = 0.99, 
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
                                     label_above_lfc,label_above_quantile,
                                     subsample_below_quantile,subsample_rate)
    dat$beta <- dat$beta.truncated
    return(ggplot_call(dat,y.label,k,max.overlaps))
  } else {

    # Create a volcano plot for each selected topic, and combine them
    # using plot_grid. This is done by recursively calling volcano_plot.
    m     <- length(k)
    plots <- vector("list",m)
    names(plots) <- k
    for (i in 1:m)
      plots[[i]] <- volcano_plot(diff_count_result,k[i],labels,y,betamax,
                                 label_above_lfc,label_above_quantile,
                                 subsample_below_quantile,subsample_rate,
                                 max.overlaps,ggplot_call,NULL)
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
                            y = c("zscore", "pvalue"), betamax = 15,
                            subsample_below_quantile = 0,
                            subsample_rate = 0.1, width = 500, height = 500,
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
                                   -Inf,0,subsample_below_quantile,
                                   subsample_rate)

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
    scale_x_continuous(expand = expansion(mult = 0.2),breaks = seq(-15,15,5)) +
    scale_y_continuous(trans = "sqrt",
      breaks = c(0,1,2,5,10,20,50,100,200,500,1e3,2e3,5e3,1e4,2e4,5e4)) +
    scale_fill_gradient2(low = "deepskyblue",mid = "gold",high = "orangered",
                         na.value = "gainsboro",
                         midpoint = mean(range(dat$mean))) +
    geom_text_repel(color = "darkgray",size = 2.25,fontface = "italic",
                    segment.color = "darkgray",segment.size = 0.25,
                    min.segment.length = 0,max.overlaps = max.overlaps,
                    na.rm = TRUE) +
    labs(x = "log-fold change",y = y.label,fill = "log10 mean",
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
  p <- plot_ly(data = dat,x = ~beta.truncated,y = ~sqrt(y),color = ~mean,
               colors = c("deepskyblue","gold","orangered"),
               text = ~sprintf(paste0("%s\nmean: %0.3f\nlogFC: %+0.3f\n",
                                      "s.e.: %0.2e\nz: %+0.3f\n",
                                      "-log10p: %0.2f"),
                               label,10^mean,beta,se,z,pval),
               type = "scatter",mode = "markers",hoverinfo = "text",
               width = width,height = height,
               marker = list(line = list(color = "white",width = 1),
                             size = 7.5))
  p <- hide_colorbar(p)
  p <- layout(p,
              xaxis = list(title = "log-fold change",zeroline = FALSE,
                           showgrid = FALSE),
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
compile_volcano_plot_data <-
  function (diff_count_result, k, labels, y, betamax, label_above_lfc,
            label_above_quantile, subsample_below_quantile, subsample_rate) {
  dat <- with(diff_count_result,
                data.frame(label = labels,
                           mean  = colmeans,
                           beta  = beta[,k],
                           se    = se[,k],
                           z     = Z[,k],
                           pval  = pval[,k],
                           y     = 0,
                           stringsAsFactors = FALSE))
  n <- nrow(dat)
  if (y == "zscore")
    dat$y <- abs(diff_count_result$Z[,k])
  else if (y == "pvalue")
    dat$y <- diff_count_result$pval[,k]
  rows     <- which(dat$mean > 0)
  dat      <- dat[rows,]
  dat$mean <- log10(dat$mean)
  if (is.infinite(label_above_quantile))
    y0 <- Inf
  else
    y0 <- quantile(dat$y,label_above_quantile)
  dat$label[with(dat,beta < label_above_lfc | y < y0)] <- ""
  if (subsample_below_quantile > 0) {
    y0    <- quantile(dat$y,subsample_below_quantile)
    rows1 <- which(dat$y >= y0)
    rows2 <- which(dat$y < y0)
    rows2 <- sample(rows2,ceiling(subsample_rate * length(rows2)))
    rows  <- sort(c(rows1,rows2))
    message(sprintf("%d out of %d data points will be included in plot",
                    length(rows),n))
    dat   <- dat[rows,]
  }
  dat <- transform(dat,beta.truncated = sign(beta) * pmin(betamax,abs(beta)))
  return(dat)
}
