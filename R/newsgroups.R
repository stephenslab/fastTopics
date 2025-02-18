#' @name newsgroups
#'
#' @title Topic modeling results from the \dQuote{20 Newsgroups} data
#'   set.
#'
#' @docType data
#' 
#' @description These are topic modeling results from the \dQuote{20
#' Newsgroups} data, with k = 10 topics. The data were originally
#' downloaded from \url{http://qwone.com/~jason/20Newsgroups} and
#' prepared by running code that found in an R Markdown file in this
#' GitHub repository:
#' \url{https://github.com/stephenslab/fastTopics-experiments}. See
#' the \dQuote{inst} directory of this package for the scripts used to
#' generate these results.
#'
#' @format \code{newsgroups} is a list with the following elements:
#' 
#' \describe{
#'
#'   \item{topics}{Original labeling of the documents: each document
#'     is from one of 20 \dQuote{newsgroups}.}
#'
#'   \item{L}{Estimated topic proportions matrix; rows are
#'     documents and columns are topics.}
#'
#'   \item{F}{Matrix containing posterior mean estimates of log-fold
#'   changes (in base-2 logarithm). These were computed using
#'   \code{\link{de_analysis}} with \code{lfc.stat = "vsnull"}. Columns
#'   are words and columns are topics.}}
#' 
#' @keywords data
#'
#' @examples
#' data(newsgroups)
#' table(newsgroups$topics)
#' dim(newsgroups$L)
#' dim(newsgroups$F)
#' 
NULL
