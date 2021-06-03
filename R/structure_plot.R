#' @title Structure Plot
#'
#' @description Create a \dQuote{Structure plot} from a multinomial topic
#'   model fit. The Structure plot represents the estimated mixture
#'   proportions of each sample as a stacked barchart, with bars of
#'   different colors representing different topics. Consequently,
#'   samples that have similar mixture proportions have similar amounts
#'   of each color.
#'
#' @details The name \dQuote{Structure plot} comes from its widespread
#' use in population genetics to visualize the results of the
#' Structure software (Rosenberg \emph{et al}, 2002).
#' 
#' For most uses of the Structure plot in population genetics, there
#' is usually some grouping of the samples (e.g., assignment to
#' pre-defined populations) that guides arrangement of the samples
#' along the horizontal axis in the bar chart. In other applications,
#' such as analysis of gene expression data, no pre-defined grouping
#' exists. Therefore, a \dQuote{smart} arrangement of the samples is,
#' by default, generated automatically by performing a 1-d t-SNE
#' embedding of the samples.
#'
#' Alternatively, a categorical variable---the grouping---may be
#' provided, in which case the samples are arranged according to that
#' grouping, then arranged within each group using t-SNE.
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
#'   \dQuote{mass} are at shown at the bottom of the plot.
#' 
#' @param colors Colors used to draw topics in Structure plot:
#'   \code{colors[1]} is the colour used to draw \code{topics[1]},
#'   \code{colors[2]} is the colour used to draw \code{topics[2]}, and
#'   so on. The default colour setting is the from
#'   \url{https://colorbrewer2.org} (qualitative data, \dQuote{9-class
#'   Set1}).
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
#'   \dQuote{sample}, \dQuote{topic} and \dQuote{mixprop}: column
#'   \dQuote{sample} contains the positions of the samples (rows of the
#'   loadings matrix) along the horizontal axis; column \dQuote{topic} is a
#'   topic (corresponding to columns of the loadings matrix); and column
#'   \dQuote{mixprop} is the mixture proportion for the given sample.
#'
#' @param ticks The placement of the group labels along the horizontal
#'   axis, and their names. For data that is not grouped, use
#'   \code{ticks = NULL}.
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
