#' @title Annotation Heatmap
#'
#' @description This is a generic plotting utility (not specific to 
#' topic the model) for comparing \dQuote{effects} across multiple
#' dimensions (e.g., topics). The function has several options for
#' selecting the features to compare.
#'
#' @param effects_matrix n x d numeric matrix, where n is the number
#'   of features and d is the number of dimensions. This could be for
#'   example the word frequencies matrix \code{F} from a multinomial
#'   topic model fitted using \code{\link{fit_topic_model}}. The row
#'   and columns of this matrix must be named, otherwise this function
#'   will throw and error.
#'
#' @param select_features This may be a character vector specifying
#'   the features to plot (rows of the effects matrix). Or it may be one
#'   of the following: \code{"largest"}, which automatically selects the
#'   largest effects for each chosen dimension; \code{"distinctive"},
#'   which automatically selects the \dQuote{most distinctive} effects
#'   for each chosen dimension; or \code{"both"}, which uses both
#'   criteria to select features. Distinctive features are defined as
#'   rows of the effects matrix that are much larger in magnitude than
#'   the effects in the other dimensions that also share the same sign.
#'
#' @param feature_sign For automated selection of features, this
#'   option determines whether to consider positive effects only
#'   (\code{"positive"}), negative effects only (\code{"negative"}), or
#'   both (\code{"both"}).
#'
#' @param dims The dimensions (columns of the effect matrix) to
#' consider for automatic feature selection. This should be dimension
#' names (not numbers).
#'
#' @param compare_dims This should be dimension names (not numbers).
#'
#' @param n For automated feature selection, the number of features to
#'   select of each type and for each dimension. (see arguments
#'   \code{select_features} and \code{feature_sign}).
#'
#' @param show_dims The dimensions (columns) to include in the
#'   plot. This should be dimension names (not numbers).
#' 
#' @param zero_value Numbers smaller than \code{zero_value} (in
#'   magnitude) are not shown in the plot.
#'
#' @param font_size Specifies the font size for the plot.
#'
#' @param verbose When \code{verbose = TRUE}, the list of selected
#'   features (rows) is printed.
#' 
#' @return A \code{ggplot} object.
#'
#' @examples
#' data(newsgroups)
#' p1 <- annotation_heatmap(newsgroups$F,feature_sign = "positive")
#' 
#' @export
#'
annotation_heatmap <-
  function (effects_matrix,
            select_features = c("largest","distinctive","both","all"),
            feature_sign = c("both","positive","negative"),
            dims = colnames(effects_matrix),
            compare_dims = colnames(effects_matrix), n = 2,
            show_dims = colnames(effects_matrix), zero_value = 0.01,
            font_size = 10, verbose = TRUE) {

  # Verify and process the effects_matrix input.
  if (!(is.matrix(effects_matrix) & is.numeric(effects_matrix)))
    stop("Input \"effects_matrix\" should be a numeric matrix")
  if (nrow(effects_matrix) < 2 | ncol(effects_matrix) < 2)
    stop("Input \"effects_matix\" should have at least 2 rows and at least ",
         "2 columns")
  if (is.null(rownames(effects_matrix)) | is.null(colnames(effects_matrix)))
    stop("Input \"effects_matrix\" should be a named matrix; ",
         "that is, the rows and columns must be named. ",
         "See help(rownames) for details.")

  # Verify and process the other inputs.
  if (!is.character(select_features))
    stop("Input \"select_features\" should be one of \"both\", ",
         "\"distinctive\", \"largest\" or \"all\", or a character ",
         "vector selecting the rows to plot")
  if (missing(select_features) |
      (length(select_features) == 1 &
       is.element(select_features[1],
                  c("both","distinctive","largest","all")))) {
    select_features <- match.arg(select_features)
    features <- NULL
  }
  else
    features <- select_features
  feature_sign <- match.arg(feature_sign)
  if (!is.character(dims))
    stop("Input \"dims\" should be a character vector specifying the ",
         "columns of the effects matrix used to select features to plot")
  if (!is.character(compare_dims))
    stop("Input \"compare_dims\" should be a character vector specifying ",
         "the columns of the effects matrix to compare to")
  if (!is.character(show_dims))
    stop("Input \"show_dims\" should be a character vector specifying the ",
         "columns of the effects matrix to include in the plot")

  # Automatically select the features to plot.
  if (is.null(features)) {
    all_features <- rownames(effects_matrix)
    if (select_features == "all")
      features <- all_features
    else {

      # Repeat for each dimension that was chosen.
      for (k in dims) {
         if (select_features == "largest")
           features <- c(features,
                         get_largest_features(effects_matrix,k,n,feature_sign))
         else if (select_features == "distinctive")
           features <-
             c(features,
               get_distinctive_features(effects_matrix,k,n,feature_sign,
                                        compare_dims))
         else if (select_features == "both")
           features <-
             c(features,
               get_largest_features(effects_matrix,k,n,feature_sign),
               get_distinctive_features(effects_matrix,k,n,feature_sign,
                                        compare_dims))
      }
    }
  }

  # Remove any selected features that were selected more than once.
  features <- features[!duplicated(features)]

  # Print the list of selected features, if requested.
  if (verbose) {
    cat("# Features selected for plot:",
        paste(features,collapse = " "),"\n")
    dput(features)
  }
  
  # Create the heatmap.
  return(effect_heatmap(effects_matrix[features,,drop = FALSE],
                        zero_value,font_size))
}

# This is the function used by annotation_heatmap() to select the
# largest positive and/or negative features for a given dimension k.
get_largest_features <-
  function (effects_matrix, k, n,
            feature_sign = c("both","positive","negative")) {
  feature_sign <- match.arg(feature_sign)
  features     <- NULL
  rows <- which(effects_matrix[,k] < 0)
  X0   <- effects_matrix[rows,,drop = FALSE]
  rows <- which(effects_matrix[,k] > 0)
  X1   <- effects_matrix[rows,,drop = FALSE]
  n0   <- min(n,nrow(X0))
  n1   <- min(n,nrow(X1))
  if ((feature_sign == "positive" | feature_sign == "both") & n1 > 0) {
    all_features <- rownames(X1)
    rows         <- order(X1[,k],decreasing = TRUE)
    rows         <- rows[1:n1]
    features     <- c(features,all_features[rows])
  }
  if ((feature_sign == "negative" | feature_sign == "both") & n0 > 0) {
    all_features <- rownames(X0)
    rows         <- order(X0[,k],decreasing = FALSE)
    rows         <- rows[1:n0]
    features     <- c(features,all_features[rows])
  }
  return(features)
}

# This is the function used by annotation_heatmap() to select the most
# distinctive positive and/or negative features for a given dimension
# k.
get_distinctive_features <-
  function (effects_matrix, k, n,
            feature_sign = c("both","positive","negative"),
            compare_dims = colnames(effects_matrix)) {
  feature_sign <- match.arg(feature_sign)
  features     <- NULL
  rows <- which(effects_matrix[,k] < 0)
  X0   <- effects_matrix[rows,,drop = FALSE]
  rows <- which(effects_matrix[,k] > 0)
  X1   <- effects_matrix[rows,,drop = FALSE]
  n0   <- min(n,nrow(X0))
  n1   <- min(n,nrow(X1))
  compare_dims <- setdiff(compare_dims,k)
  if ((feature_sign == "positive" | feature_sign == "both") & n1 > 0) {
    all_features <- rownames(X1)
    y            <- X1[,k] - apply(X1[,compare_dims,drop = FALSE],1,max)
    rows         <- order(y,decreasing = TRUE)
    rows         <- rows[1:n1]
    features     <- c(features,all_features[rows])
  }
  if ((feature_sign == "negative" | feature_sign == "both") & n0 > 0) {
    all_features <- rownames(X0)
    y            <- X0[,k] - apply(X0[,compare_dims],1,min)
    rows         <- order(y,decreasing = FALSE)
    rows         <- rows[1:n0]
    features     <- c(features,all_features[rows])
  }
  return(features)
}

# This function is used by annotation_heatmap() to create the final heatmap.
#
#' @importFrom stats quantile
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 scale_size
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggplot2 guides
#' @importFrom ggplot2 guide_legend
#' @importFrom ggplot2 labs
#' @importFrom cowplot theme_cowplot
#' 
effect_heatmap <- function (effects_matrix, zero_value, font_size) {
  features <- rownames(effects_matrix)
  pdat <- data.frame(feature_name = features,stringsAsFactors = FALSE)
  pdat <- cbind(pdat,effects_matrix)
  pdat <- melt(pdat,id.vars = "feature_name",variable.name = "dim",
               value.name = "value")
  pdat <- transform(pdat,
                    effect_size = abs(value),
                    effect_sign = factor(sign(value),c(-1,0,1)),
                    feature_name = factor(feature_name,rev(features)))
  pdat$effect_size[abs(pdat$effect_size) < zero_value] <- NA
  effect_size_breaks <-
    unname(quantile(pdat$effect_size,probs = c(0,0.25,0.5,0.75,1),
                    na.rm = TRUE))
  if (any(pdat$effect_sign == -1))
    dot_colors <- c("navy","gray","orangered")
  else
    dot_colors <- c("navy","gray","black")
  return(ggplot(pdat,aes(x = dim,y = feature_name,size = effect_size,
                         fill = effect_sign)) +
         geom_point(color = "white",shape = 21,na.rm = TRUE) +
         scale_size(range = c(1,6),breaks = effect_size_breaks,
                    labels = round(effect_size_breaks,digits = 3)) +
         scale_fill_manual(values = dot_colors,drop = FALSE) +
         guides(size = guide_legend(override.aes = list(fill = "black")),
                fill = guide_legend(override.aes = list(color = "white",
                                                        shape = 21,
                                                        size = 3))) +
         labs(x = "",y = "",fill = "effect sign",size = "effect size") +
         theme_cowplot(font_size = font_size))
}
