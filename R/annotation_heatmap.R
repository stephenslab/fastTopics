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
            select_features = c("both", "distinctive", "largest"),
            feature_sign = c("both", "positive", "negative"),
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
    stop("Input \"select_features\" should be one of \"both\"",
         "\"distinctive\" or \"largest\", or a character vector selecting
         the rows to plot")
  if (length(select_features) == 1 &
      is.element(select_features,c("both","distinctive","largest")))
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
