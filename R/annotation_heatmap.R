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
#' @return TO DO.
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
            compare_dims = dims) {

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
      is.element(select_features,c("both","distinctive","largest","all")))
    select_features <- match.arg(select_features)
  feature_sign <- match.arg(feature_sign)
  if (!is.character(dims))
    stop("Input \"dims\" should be a character vector specifying the ",
         "columns of the effects matrix to plot")
  if (!is.character(compare_dims))
    stop("Input \"compare_dims\" should be a character vector specifying ",
         "the columns of the effects matrix to compare to")

  
  
  return(NULL)
}

# TO DO: Explain here what this function is for, and how it is used.
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
effect_plot <- function (effects_matrix, zero_value = 0.01, font_size = 10) {
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
