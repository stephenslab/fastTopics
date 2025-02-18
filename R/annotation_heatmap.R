#' @title Annotation Heatmap
#'
#' @description TO DO.
#'
#' @param effects_matrix Describe the effects_matrix input argument here.
#'
#' @param select_features Describe the select_features input argument here.
#'
#' @param feature_sign Describe the feature_sign input argument here.
#'
#' @param dims Describe the dims input argument here.
#'
#' @param compare_dims Describe the compare_dims argument here.
#'
#' @param n Describe argument n here.
#'
#' @param zero_value Describe the zero_value argument here.
#'
#' @param font_size Describe the font_size argument here.
#' 
#' @return A \code{ggplot} object.
#'
#' @examples
#' # TO DO.
#' 
#' @export
#'
annotation_heatmap <-
  function (effects_matrix,
            select_features = c("both","distinctive","largest","all"),
            feature_sign = c("both","positive","negative"),
            dims = colnames(effects_matrix),
            compare_dims = dims,
            n = 2,
            zero_value = 0.01,
            font_size = 10) {

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
  if (length(select_features) == 1 &
      is.element(select_features[1],c("both","distinctive","largest","all"))) {
    select_features <- match.arg(select_features)
    features <- NULL
  }
  else
    features <- select_features
  feature_sign <- match.arg(feature_sign)
  if (!is.character(dims))
    stop("Input \"dims\" should be a character vector specifying the ",
         "columns of the effects matrix to plot")
  if (!is.character(compare_dims))
    stop("Input \"compare_dims\" should be a character vector specifying ",
         "the columns of the effects matrix to compare to")

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
         else if (select_features == "distinctive") {
           # Add code here.
         } else if (select_features == "both")
           features <- c(features,
                         get_largest_features(effects_matrix,k,n,feature_sign))
      }
    }
  }

  # Remove any selected features that were selected more than once.
  features <- features[!duplicated(features)]
  
  # Create the heatmap.
  return(effect_heatmap(effects_matrix[features,dims],zero_value,font_size))
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
  pdat$effect_size[pdat$effect_size < zero_value] <- NA
  effect_size_breaks <- quantile(pdat$effect_size,probs = c(0,0.25,0.5,0.75,1),
                                 na.rm = TRUE)
  effect_size_breaks <- unname(round(effect_size_breaks,digits = 2))
  return(ggplot(pdat,aes(x = dim,y = feature_name,size = effect_size,
                         fill = effect_sign)) +
         geom_point(color = "white",shape = 21,na.rm = TRUE) +
         scale_size(range = c(1,6),breaks = effect_size_breaks) +
         scale_fill_manual(values = c("navy","gray","orangered"),
                           drop = FALSE) +
         guides(size = guide_legend(override.aes = list(fill = "black")),
                fill = guide_legend(override.aes = list(color = "white",
                                                        shape = 21,
                                                        size = 3))) +
         labs(x = "",y = "",fill = "effect sign",size = "effect size") +
         theme_cowplot(font_size = font_size))
}
