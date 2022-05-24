#' @title Combine Topics in Multinomial Topic Model
#'
#' @description Combine two or more topics in a multinomial topic
#'   model fit.
#'
#' @details Mixture proportions are combined by summation, and factors
#'   are combined by averaging.
#' 
#' @param fit A multinomial topic model fit.
#'
#' @param k The names or numbers of the topics to be combined. Two or
#'   more topics should be chosen.
#'
#' @return A multinomial topic model fit.
#' 
#' @export
#' 
merge_topics <- function (fit, k) {

  # Verify input "fit".
  if (!inherits(fit,"multinom_topic_model_fit"))
    stop("Input argument \"fit\" should be an object of class ",
         "\"multinom_topic_model_fit\"")
  verify.fit(fit)
  
  # Verify and process input "k".
  msg <- paste("Input argument \"k\" should contain valid topic names or",
               "numbers (column indices of F and L)")
  if (!((is.numeric(k) | is.character(k)) & length(k) >= 2))
    stop(msg)
  if (is.numeric(k)) {
    if (!all(k >= 1 & k <= ncol(fit$F)))
      stop(msg)
  } else {
    if (!all(is.element(k,colnames(fit$F))))
      stop(msg)
    k <- match(k,colnames(fit$F))
  }
  
  # Combine the selected topics.
  out1   <- combine_factors(fit$F,fit$L,k)
  out2   <- combine_factors(fit$Fn,fit$Ln,k)
  out3   <- combine_factors(fit$Fy,fit$Ly,k)
  fit$F  <- out1$F
  fit$L  <- out1$L
  fit$Fn <- out2$F
  fit$Ln <- out2$L
  fit$Fy <- out3$F
  fit$Ly <- out3$L
  return(fit)
}

# Combine two or more columns of the factors matrix (F) and loadings
# matrix (L). Loadings are combined by summation, and factors are
# combined by averaging.
combine_factors <- function (F, L, k) {
  if (is.null(colnames(F)))
    y <- NULL
  else {
    y <- colnames(F)
    y <- c(y[-k],paste(y[k],collapse = "+"))
  }
  F <- cbind(F[,-k],rowMeans(F[,k]))
  L <- cbind(L[,-k],rowSums(L[,k]))
  colnames(F) <- y
  colnames(L) <- y
  return(list(F = F,L = L))
}
