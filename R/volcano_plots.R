#' @rdname volcano_plot
#'
#' @title Volcano Plot
#'
#' @description Create a \dQuote{volcano} plot to visualize the
#'   results of a differential count analysis using a topic model. Here,
#'   the volcano plot is a scatterplot in which the posterior mean
#'   log-fold change (LFC), estimated by running the methods implemented
#'   in \code{\link{de_analysis}}, is plotted against the estimated
#'   z-score.
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
#' package. The \dQuote{hover text} shows the label (see input
#' argument \dQuote{labels}) and detailed LFC statistics as they were
#' calculated by \code{\link{de_analysis}}.
#' 
#' @param de An object of class \dQuote{topic_model_de_analysis},
#' usually an output from \code{\link{de_analysis}}.
#'
#' @param k The topic, selected by number or name.
#'
#' @param labels Character vector specifying how the points in the
#'   volcano plot are labeled. This should be a character vector with
#'   one entry per LFC estimate (row of \code{de$postmean}). When not
#'   specified, the row names of \code{de$postmean} are used. When
#'   available. labels are added to the plot using
#'   \code{\link[ggrepel]{geom_text_repel}}.
#'
#' @param do.label Describe input argument "do.label" here.
#' 
#' @param ymax Describe input argument "ymax" here.
#' 
#' @param max.overlaps Argument passed to
#'   \code{\link[ggrepel]{geom_text_repel}}.
#' 
#' @param ggplot_call The function used to create the plot. Replace
#'   \code{volcano_plot_ggplot_call} with your own function to customize
#'   the appearance of the plot.
#'
#' @return A \code{ggplot} object or a \code{plotly} object.
#'
#' @seealso \code{\link{de_analysis}}
#' 
#' @examples
#' # See help(de_analysis) for examples.
#' 
#' @export
#' 
volcano_plot <-
  function (de, k, labels, do.label = volcano_plot_do_label_default,
            ymax = Inf, max.overlaps = Inf,
            ggplot_call = volcano_plot_ggplot_call) {
  if (!inherits(de,"topic_model_de_analysis"))
    stop("Input \"de\" should be an object of class ",
         "\"topic_model_de_analysis\"")
  n <- nrow(de$postmean)
  if (missing(labels)) {
    if (!is.null(rownames(de$postmean)))
      labels <- rownames(de$postmean)
    else
      labels <- as.character(seq(1,n))
  }
  if (!(is.character(labels) & length(labels) == n))
    stop("Input argument \"labels\", when specified, should be a character ",
         "vector with one entry per log-fold change estimate (column of ",
         "the counts matrix)")
  dat <- compile_volcano_plot_data(de,k,labels,do.label)
  return(ggplot_call(dat,max.overlaps))
}

#' @rdname volcano_plot
#'
#' @param x An object of class \dQuote{topic_model_de_analysis},
#'   usually an output from \code{\link{de_analysis}}.
#'
#' @param \dots Additional arguments passed to \code{volcano_plot}.
#' 
#' @importFrom graphics plot
#' 
#' @method plot topic_model_de_analysis
#'
#' @export
#'
plot.topic_model_de_analysis <- function (x, ...)
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
volcano_plotly <- function (de, k, file, labels,
                            subsample_below_quantile = 0,
                            subsample_rate = 0.1, width = 500, height = 500,
                            title = paste("topic",k),
                            plot_ly_call = volcano_plot_ly_call) {

  # Check and process input arguments.
  y <- match.arg(y)
  if (!inherits(de_result,"topic_model_de_analysis"))
    stop("Input \"de\" should be an object of class ",
         "\"topic_model_de_analysis\"")
  beta <- de_result$beta
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
  dat <- compile_volcano_plot_data(de_result,k,labels,y,betamax,
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
#' @param lfc Describe input argument "lfc" here.
#'
#' @param z Describe input argument "z" here.
#'
#' @importFrom stats quantile
#' 
#' @export
#' 
volcano_plot_do_label_default <- function (lfc, z)
  abs(z) >= quantile(abs(z),0.999,na.rm = TRUE) |
    lfc <= quantile(lfc,0.001,na.rm = TRUE) |
    lfc >= quantile(lfc,0.999,na.rm = TRUE)

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
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_text
#' @importFrom ggrepel geom_text_repel
#' @importFrom cowplot theme_cowplot
#' 
#' @export
#' 
volcano_plot_ggplot_call <- function (dat, max.overlaps = Inf, font.size = 9)
  ggplot(dat,aes_string(x="postmean",y = "z",fill = "lfsr",label = "label")) +
    geom_point(color = "white",stroke = 0.3,shape = 21,na.rm = TRUE) +
    geom_text_repel(color = "darkgray",size = 2.25,fontface = "italic",
                    segment.color = "darkgray",segment.size = 0.25,
                    min.segment.length = 0,max.overlaps = max.overlaps,
                    na.rm = TRUE) +
    #scale_x_continuous(expand=expansion(mult = 0.2),breaks = seq(-15,15,5)) +
    scale_y_continuous(trans = "sqrt",
      breaks = c(0,1,2,5,10,20,50,100,200,500,1e3,2e3,5e3,1e4,2e4,5e4)) +
    # scale_fill_gradient2(low = "deepskyblue",mid = "gold",high = "orangered",
    #                      na.value = "gainsboro",
    #                      midpoint = mean(range(dat$mean))) +
    labs(x = "log-fold change",y = "|z-score|") +
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
compile_volcano_plot_data <- function (de, k, labels, do.label) {
  if (is.null(de$lfsr))
    lfsr <- as.numeric(NA)
  else
    lfsr <- de$lfsr[,k]
  dat <- data.frame(label    = labels,
                    postmean = de$postmean[,k],
                    z        = de$z[,k],
                    lfsr     = lfsr,
                    stringsAsFactors = FALSE)
  dat$label[which(!do.label(dat$postmean,dat$z))] <- ""
  dat$z <- abs(dat$z)
  return(dat)
}
