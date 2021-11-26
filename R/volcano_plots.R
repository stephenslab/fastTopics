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
#' @details Interactive volcano plots can be created using the
#'   \dQuote{plotly} package. The \dQuote{hover text} shows the label
#'   (see input argument \dQuote{labels}) and detailed LFC statistics as
#'   they were calculated by \code{\link{de_analysis}}.
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
#' @param do.label The function used to deetermine which LFC estimates
#'   to label. Replace \code{volcano_plot_do_label_default} with your
#'   own function to customize the labeling of points in the volcano
#'   plot.
#' 
#' @param ymax z-scores greater than \code{ymax} (in magnitude) are
#'   shown as \code{max}. Setting \code{ymax} to a finite value can
#'   improve the volcano plot when some z-scores are much larger (in
#'   magnitude) than others.
#' 
#' @param max.overlaps Argument passed to
#'   \code{\link[ggrepel]{geom_text_repel}}.
#'
#' @param plot.title The title of the plot.
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
            ymax = Inf, max.overlaps = Inf, plot.title = paste("topic",k),
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
  dat <- compile_volcano_plot_data(de,k,ymax,labels,do.label)
  return(ggplot_call(dat,plot.title,max.overlaps))
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
#' @param lfc A vector of log-fold change estimates.
#'
#' @param z A vector of z-scores of the same length as \code{lfc}.
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
#'   \dQuote{postmean}, \dQuote{z}, \dQuote{lfsr} and \dQuote{label}.
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
volcano_plot_ggplot_call <- function (dat, plot.title, max.overlaps = Inf,
                                      font.size = 9) {
  dat$lfsr <- cut(dat$lfsr,c(-1,0.001,0.01,0.05,Inf))
  return(ggplot(dat,aes_string(x = "postmean",y = "z",fill = "lfsr",
                               label = "label")) +
         geom_point(color = "white",stroke = 0.3,shape = 21,na.rm = TRUE) +
         geom_text_repel(color = "darkgray",size = 2.25,fontface = "italic",
                         segment.color = "darkgray",segment.size = 0.25,
                         min.segment.length = 0,max.overlaps = max.overlaps,
                         na.rm = TRUE) +
         scale_y_continuous(trans = "sqrt",
                            breaks = c(0,1,2,5,10,20,50,100,200,500,
                                       1000,2000,5000,1e4,2e4,5e4)) +
         scale_fill_manual(values = c("deepskyblue","gold","orange","coral")) +
         labs(x = "log-fold change",y = "|z-score|",title = plot.title) +
         theme_cowplot(font.size) +
         theme(plot.title = element_text(size = font.size,face = "plain")))
}

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
compile_volcano_plot_data <- function (de, k, ymax, labels, do.label) {
  dat <- data.frame(label    = labels,
                    postmean = de$postmean[,k],
                    z        = de$z[,k],
                    lfsr     = de$lfsr[,k],
                    stringsAsFactors = FALSE)
  dat$label[which(!do.label(dat$postmean,dat$z))] <- ""
  dat$z <- pmin(ymax,abs(dat$z))
  return(dat)
}
