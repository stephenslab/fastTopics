#' @name droplet
#'
#' @title Droplet single-cell RNA-seq read count data from Montoro
#'   \emph{et al} (2018)
#'
#' @docType data
#' 
#' @description These data are gene expression profiles of trachea
#' epithelial cells in C57BL/6 mice obtained using droplet-based 3'
#' single-cell RNA-seq. They were prepared from file
#' \code{GSE103354_Trachea_droplet_UMIcounts.txt.gz} downloaded from
#' the Gene Expression Omnibus (GEO) website, accession GSE103354.
#' 
#' @format \code{droplet} is a 7,193 x 17,133 sparse matrix of read
#' counts, with rows corresponding to samples (cells), and columns
#' corresponding to genes.
#' 
#' @references
#'
#' D. T. Montoro \emph{et al} (2018). A revised airway epithelial
#' hierarchy includes CFTR-expressing ionocytes. \emph{Nature} \bold{560},
#' 319–-324.
#' 
#' @keywords data
#'
#' @examples
#'
#' # Roughly 10% of the read counts are greater than zero.
#' data(droplet)
#' nnzero(droplet)/length(droplet)
#' 
NULL
