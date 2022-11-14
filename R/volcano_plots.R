#' @rdname volcano_plot
#'
#' @title Volcano Plots for Visualizing Results of Differential Expression Analysis
#'
#' @description Create a \dQuote{volcano} plot to visualize the
#'   results of a differential count analysis using a topic model. Here,
#'   the volcano plot is a scatterplot in which the posterior mean
#'   log-fold change (LFC), estimated by running the methods implemented
#'   in \code{\link{de_analysis}}, is plotted against the estimated
#'   z-score. Use \code{volcano_plotly} to create an interactive volcano
#'   plot.
#'
#' @details Interactive volcano plots can be created using the
#'   \sQuote{plotly} package. The \dQuote{hover text} shows the label
#'   and detailed LFC statistics.
#' 
#' @param de An object of class \dQuote{topic_model_de_analysis},
#'   usually an output from \code{\link{de_analysis}}. It is better to
#'   run \code{de_analysis} with \code{shrink.method = "ash"} so that
#'   the points in the volcano plot can be coloured by their local false
#'   sign rate (lfsr).
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
#' @param plot_ly_call The function used to create the plot. Replace
#'   \code{volcano_plot_ly_call} with your own function to customize
#'   the appearance of the interactive plot.
#'
#' @importFrom htmlwidgets saveWidget
#' 
#' @export
#' 
volcano_plotly <- function (de, k, file, labels,
                            ymax = Inf, width = 500, height = 500,
                            plot.title = paste("topic",k),
                            plot_ly_call = volcano_plot_ly_call) {

  # Check and process input arguments.
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
  dat <- compile_volcano_plot_data(de,k,ymax,labels)
  p <- volcano_plot_ly_call(dat,plot.title,width,height)
  if (!missing(file))
    saveWidget(p,file,selfcontained = TRUE,title = plot.title)
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
#'   \dQuote{f0}, \dQuote{postmean}, \dQuote{z}, \dQuote{zabs},
#'   \dQuote{lfsr} and \dQuote{label}.
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
  return(ggplot(dat,aes_string(x = "postmean",y = "zabs",fill = "lfsr",
                               label = "label")) +
         geom_point(color = "white",stroke = 0.3,shape = 21,na.rm = TRUE) +
         geom_text_repel(color = "darkgray",size = 2.25,fontface = "italic",
                         segment.color = "darkgray",segment.size = 0.25,
                         min.segment.length = 0,max.overlaps = max.overlaps,
                         na.rm = TRUE) +
         scale_y_continuous(trans = "sqrt",
                            breaks = c(0,1,2,5,10,20,50,100,200,500,
                                       1000,2000,5000,1e4,2e4,5e4)) +
         # scale_y_continuous(trans = "log10") +
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
volcano_plot_ly_call <- function (dat, plot.title, width, height) {
  dat$fill <- cut(dat$lfsr,c(-1,0.001,0.01,0.05,Inf))
  p <- plot_ly(data = dat,x = ~postmean,y = ~sqrt(zabs),color = ~fill,
               colors = c("deepskyblue","gold","orange","coral"),
               text = ~sprintf(paste0("%s\nmean(null): %0.2e\n",
                                      "logFC(lower): %+0.3f\n",
                                      "logFC(mean): %+0.3f\n",
                                      "logFC(upper): %0.3f\n",
                                      "z-score: %+0.3f\nlfsr: %0.2e"),
                               label,f0,lower,postmean,upper,z,lfsr),
               type = "scatter",mode = "markers",hoverinfo = "text",
               width = width,height = height,
               marker = list(line = list(color = "white",width = 1),
                             size = 7.5))
  p <- hide_colorbar(p)
  p <- layout(p,
              xaxis = list(title = "log-fold change",zeroline = FALSE,
                           showgrid = FALSE),
              yaxis = list(title = "sqrt|z-score|",zeroline = FALSE,
                           showgrid = FALSE),
              hoverlabel = list(bgcolor = "white",bordercolor = "black",
                                font = list(color = "black",family = "arial",
                                            size = 12)),
              font = list(family = "arial",size = 12),
              showlegend = FALSE,title = plot.title)
  return(p)
}

# This is used by volcano_plot and volcano_plotly to compile the data
# frame passed to 'ggplot'.
compile_volcano_plot_data <- function (de, k, ymax, labels, do.label = NULL) {
  if (all(is.na(de$lfsr))) {
    message("lfsr is not available, probably because \"shrink.method\" was ",
            "not set to \"ash\"; lfsr in plot should be ignored")
    lfsr <- 0
  } else
    lfsr <- de$lfsr[,k]
  dat <- data.frame(label    = labels,
                    f0       = de$f0,
                    lower    = de$lower[,k],
                    postmean = de$postmean[,k],
                    upper    = de$upper[,k],
                    z        = de$z[,k], # fit$F[,k] + 1e-5,
                    lfsr     = lfsr,
                    stringsAsFactors = FALSE)
  if (!is.null(do.label))
    dat$label[which(!do.label(dat$postmean,dat$z))] <- ""
  dat$zabs <- pmin(ymax,abs(dat$z))
  return(dat)
}
