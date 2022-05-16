#' @title Structure Plot
#'
#' @description Create a \dQuote{Structure plot} from a multinomial topic
#'   model fit. The Structure plot represents the estimated topic
#'   proportions of each sample in a stacked bar chart, with bars of
#'   different colors representing different topics. Consequently,
#'   samples that have similar topic proportions have similar amounts
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
#' such as analysis of gene expression data, a pre-defined grouping
#' may not always be available. Therefore, a \dQuote{smart}
#' arrangement of the samples is, by default, generated automatically
#' by performing a 1-d embedding of the samples.
#'
#' @param fit An object of class \dQuote{poisson_nmf_fit} or
#'   \dQuote{multinom_topic_model_fit}. If a Poisson NMF fit is provided
#'   as input, the corresponding multinomial topic model fit is
#'   automatically recovered using \code{\link{poisson2multinom}}.
#'
#' @param topics Top-to-bottom ordering of the topics in the Structure
#'   plot; \code{topics[1]} is shown on the top, \code{topics[2]} is
#'   shown next, and so on. If the ordering of the topics is not
#'   specified, the topics are automatically ordered so that the topics
#'   with the greatest total \dQuote{mass} are at shown at the bottom of
#'   the plot. The topics may be specified by number or by name.
#' 
#' @param grouping Optional categorical variable (a factor) with one
#'   entry for each row of the loadings matrix \code{fit$L} defining a
#'   grouping of the samples (rows). The samples (rows) are arranged
#'   along the horizontal axis according to this grouping, then within
#'   each group according to \code{loadings_order}. If
#'   \code{grouping} is not a factor, an attempt is made to convert it
#'   to a factor using \code{as.factor}. Note that if
#'   \code{loadings_order} is specified manually, \code{grouping} should
#'   be the groups for the rows of \code{fit$L} \emph{before} reordering.
#' 
#' @param loadings_order Ordering of the rows of the loadings matrix
#'   \code{fit$L} along the horizontal axis the Structure plot (after
#'   they have been grouped). If \code{loadings_order = "embed"}, the
#'   ordering is generated automatically from a 1-d embedding,
#'   separately for each group. The rows may be specified by number or
#'   by name. Note that \code{loadings_order} may include all the rows
#'   of \code{fit$L}, or a subset.
#'
#' @param n The maximum number of samples (rows of the loadings matrix
#'   \code{fit$L}) to include in the plot. Typically there is little to
#'   no benefit in including large number of samples in the Structure
#'   plot due to screen resolution limits. Ignored if
#'   \code{loadings_order} is provided.
#'
#' @param colors Colors used to draw topics in Structure plot. The
#'   default colour setting is the from \url{https://colorbrewer2.org}
#'   (qualitative data, \dQuote{9-class Set1}).
#'
#' @param gap The horizontal spacing between groups. Ignored if
#'   \code{grouping} is not provided.
#'
#' @param embed_method The function used to compute an 1-d embedding
#'   from a loadings matrix \code{fit$L}; only used if
#'   \code{loadings_order = "embed"}. The function must accept the
#'   multinomial topic model fit as its first input (\dQuote{fit}) and
#'   additional arguments may be passed (\dots). The output should be a
#'   named numeric vector with one entry per row of \code{fit$L}, and
#'   the names of the entries should be the same as the row names of
#'   \code{fit$L}.
#'
#' @param ggplot_call The function used to create the plot. Replace
#'   \code{structure_plot_ggplot_call} with your own function to
#'   customize the appearance of the plot.
#'
#' @param \dots Additional arguments passed to \code{structure_plot}
#'   (for the \code{plot} method) or \code{embed_method} (for
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
#' @examples
#' set.seed(1)
#' data(pbmc_facs)
#'
#' # Get the multinomial topic model fitted to the
#' # PBMC data.
#' fit <- pbmc_facs$fit
#'
#' # Create a Structure plot without labels. The samples (rows of L) are
#' # automatically arranged along the x-axis using t-SNE to highlight the
#' # structure in the data.
#' p1 <- structure_plot(fit)
#'
#' # Create a Structure plot with the FACS cell-type labels. Within each
#' # group (cell-type), the cells (rows of L) are automatically arranged
#' # using t-SNE.
#' subpop <- pbmc_facs$samples$subpop
#' p2 <- structure_plot(fit,grouping = subpop)
#'
#' # Next, we apply some customizations to improve the plot: (1) use the
#' # "topics" argument to specify the order in which the topic
#' # proportions are stacked on top of each other; (2) use the "gap"
#' # argrument to increase the whitespace between the groups; (3) use "n"
#' # to decrease the number of rows of L included in the Structure plot;
#' # and (4) use "colors" to change the colors used to draw the topic
#' # proportions.
#' topic_colors <- c("skyblue","forestgreen","darkmagenta",
#'                   "dodgerblue","gold","darkorange")
#' p3 <- structure_plot(fit,grouping = pbmc_facs$samples$subpop,gap = 20,
#'                      n = 1500,topics = c(5,6,1,4,2,3),colors = topic_colors)
#'
#' # In this example, we use UMAP instead of t-SNE to arrange all 3,744
#' # cells in the Structure plot. Note that this can be accomplished in a
#' # different way by overriding the default setting of "embed_method".
#' y <- drop(umap_from_topics(fit,dims = 1))
#' p4 <- structure_plot(fit,loadings_order = order(y),grouping = subpop,
#'                      gap = 40,colors = topic_colors)
#'
#' # In this final example, we plot a random subset of 400 cells, and
#' # arrange the cells randomly along the horizontal axis of the
#' # Structure plot.
#' p5 <- structure_plot(fit,loadings_order = sample(3744,400),gap = 10,
#'                      grouping = subpop,colors = topic_colors)
#'
#' @export
#'
structure_plot <-
  function (fit, topics, grouping, loadings_order = "embed", n = 2000,
            colors = c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00",
                       "#ffff33","#a65628","#f781bf","#999999"),
            gap = 1, embed_method = structure_plot_default_embed_method,
            ggplot_call = structure_plot_ggplot_call, ...) {

  # Check and process input argument "fit".
  if (!(inherits(fit,"poisson_nmf_fit") |
        inherits(fit,"multinom_topic_model_fit")))
    stop("Input \"fit\" should be an object of class \"poisson_nmf_fit\" or ",
         "\"multinom_topic_model_fit\"")
  if (inherits(fit,"poisson_nmf_fit"))
    fit <- poisson2multinom(fit)
  n0 <- nrow(fit$L)
  k  <- ncol(fit$L)
  
  # Check and process input argument "topics".
  if (is.null(colnames(fit$L)))
    colnames(fit$L) <- paste0("k",1:k)
  if (missing(topics))
    topics <- order(colMeans(fit$L))
  if (!is.character(topics))
    topics <- colnames(fit$L[,topics,drop = FALSE])
  if (!(length(topics) > 1 & all(is.element(topics,colnames(fit$L)))))
    stop("Input argument \"topics\" should be a subset of at least two ",
         "topics (columns of fit$L) specified by their names or column ",
         "indices")

  # Check and process input argument "grouping".
  if (missing(grouping))
    grouping <- factor(rep(1,n0))
  if (!is.factor(grouping))
    grouping <- as.factor(grouping)
  if (length(grouping) != n0)
    stop("Input argument \"grouping\" should be a factor with one entry ",
         "for each row of fit$L")

  # Check and process input argument "colors".
  if (length(colors) < k)
    stop("There must be at least as many colours as topics")
  names(colors) <- colnames(fit$L)
  colors <- colors[topics]

  # Check and process input arguments "loadings_order" and "n".
  if (all(loadings_order == "embed")) {

    # If necessary, randomly subsample the rows of L.
    if (n < n0) {
      rows <- sample(n0,n)
      fit <- select_loadings(fit,rows)
      grouping <- grouping[rows,drop = FALSE]
    }

    # The ordering of the rows is not provided, so determine an
    # ordering by computing a 1-d embedding of L.
    if (nlevels(grouping) == 1) {
      y <- embed_method(fit,...)
      loadings_order <- order(y)
    } else {
      loadings_order <- NULL
      for (group in levels(grouping)) {
        i <- which(grouping == group)
        if (length(i) > 0)
          y <- embed_method(select_loadings(fit,i),...)
        loadings_order <- c(loadings_order,i[order(y)])
      }
    }
  } else {
    if (!missing(n))
      warning("Input argument \"n\" is ignored when \"loadings_order\" is ",
              "not \"embed\"")
    if (is.character(loadings_order))
      loadings_order <- match(loadings_order,rownames(fit$L))
  }
  
  # Prepare the data for plotting and create the structure plot.
  fit$L <- fit$L[loadings_order,]
  grouping <- grouping[loadings_order,drop = TRUE]
  if (nlevels(grouping) == 1) {
    dat <- compile_structure_plot_data(fit$L,topics)
    return(ggplot_call(dat,colors))
  } else {
    out <- compile_grouped_structure_plot_data(fit$L,topics,grouping,gap)
    return(ggplot_call(out$dat,colors,out$ticks))
  }
}

#' @rdname structure_plot
#'
#' @importFrom stats rnorm
#' 
#' @export
#' 
structure_plot_default_embed_method <- function (fit,...) {
  if (nrow(fit$L) < 20)
    return(rnorm(nrow(fit$L)))
  else {
    d <- dim(fit$L)
    message(sprintf("Running tsne on %s x %s matrix.",d[1],d[2]))
    return(drop(suppressMessages(tsne_from_topics(fit,dims = 1,...))))
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
#'   \dQuote{sample}, \dQuote{topic} and \dQuote{prop}: the
#'   \dQuote{sample} column contains the positions of the samples (rows
#'   of the L matrix) along the horizontal axis; the \dQuote{topic}
#'   column is a topic (a column of L); and the \dQuote{prop} column is
#'   the topic proportion for the respective sample.
#'
#' @param ticks The placement of the group labels along the horizontal
#'   axis, and their names. For data that are not grouped, use
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
  ggplot(dat,aes_string(x = "sample",y = "prop",color = "topic",
                        fill = "topic")) +
    geom_col() +
    scale_x_continuous(limits = c(0,max(dat$sample) + 1),breaks = ticks,
                       labels = names(ticks)) +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors) +
    labs(x = "",y = "topic proportion") +
    theme_cowplot(font.size) +
    theme(axis.line   = element_blank(),
          axis.ticks  = element_blank(),
          axis.text.x = element_text(angle = 45,hjust = 1))

# This is used by structure_plot to create a data frame suitable for
# plotting with ggplot. Input argument L is the topic proportions
# matrix. Input argument "topics" is the vector of the selected topics
# (that is, selected columns of L). The output is a data frame with
# three columns: "sample", a row of L (numeric); "topic", a topic
# (factor); and "prop", the topic proportion for the given sample
# (numeric).
compile_structure_plot_data <- function (L, topics) {
  n <- nrow(L)
  k <- length(topics)
  dat <- data.frame(sample = rep(1:n,times = k),
                    topic  = rep(topics,each = n),
                    prop   = c(L[,topics]))
  dat$topic <- factor(dat$topic,topics)
  return(dat)
}

# This is used by structure_plot to create a data frame suitable for
# plotting with ggplot when the data are grouped. Input argument L is
# the topic proportions matrix. Input argument "topics" is the vector
# of selected topics (that is, selected columns of L).  Input argument
# "grouping" is a factor with one entry for each row of L giving the
# assigned group. The rows of L (and, correspondingly, the grouping
# vector) should already ordered by the groups; that is, grouping =
# sort(grouping). Finally, a "gap" is added to the sample indices in
# each group to provide a visual spacing of the groups.
compile_grouped_structure_plot_data <- function (L, topics, grouping,
                                                 gap = 0) {
  ticks <- rep(0,nlevels(grouping))
  names(ticks) <- levels(grouping)
  dat <- NULL
  m <- 0
  for (j in levels(grouping)) {
    i          <- which(grouping == j)
    out        <- compile_structure_plot_data(L[i,,drop = FALSE],topics)
    out$sample <- out$sample + m
    n          <- length(i)
    dat        <- rbind(dat,out)
    ticks[j]   <- m + n/2
    m          <- m + n + gap
  }
  return(list(dat = dat,ticks = ticks))
}
