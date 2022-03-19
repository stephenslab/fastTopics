#' @title Perform HOMER Motif Enrichment Analysis using DE Genomic Positions
#'
#' @description Run HOMER motif finding algorithm
#'   (\code{findMotifsGenome.pl}) to identify motifs enriched for
#'   differentially expressed (DE) genomic positions. See
#'   \url{http://homer.ucsd.edu} for more information.
#'
#' @param de An object of class \dQuote{topic_model_de_analysis},
#'   usually the result of running \code{\link{de_analysis}}.
#'
#' @param k Use the DE analysis results for this topic.
#'
#' @param positions A table of genomic positions corresponding to rows
#'   of the \code{de_analysis} results. Specifically, it should a data
#'   frame with four columns: \dQuote{chr}, chromosome name or number;
#'   \dQuote{start}, start position of genomic feature; \dQuote{end},
#'   end position of genomic feature; and \dQuote{name}, the name of the
#'   genomic feature. If not specified, the genomic positions will be
#'   extracted from the row names of \code{de$postmean}, in which the
#'   row names are expected to be of the form \code{chr_start_end}. The
#'   genomic positions will be written to a BED file (see
#'   \url{https://genome.ucsc.edu/FAQ/FAQformat.html} for more
#'   information about BED files).
#'  
#' @param genome The genome parameter passed to
#'   \code{findMotifsGenome.pl}.
#' 
#' @param subset Describe input argument "subset" here.
#'
#' @param homer.exec The name or file path of the HOMER
#'   \code{findMotifsGenome.pl} excutable.
#' 
#' @param out.dir The positions BED file and HOMER results are written
#'   to this directory.
#'
#' @param homer.options Character string used to override default
#'   \code{findMotifsGenome.pl} options.
#'
#' @param verbose When \code{verbose = TRUE}, progress information is
#'   printed to the console.
#' 
#' @return A data frame containing the motif enrichment results. It
#'   is created from the \code{knownResults.txt} HOMER output.
#'
#' @importFrom utils read.table
#' @importFrom utils write.table
#' 
#' @references
#' Heinz, S., Benner, C., Spann, N., Bertolino, E., Lin, Y. C., Laslo,
#' P., Cheng, J. X., Murre, C., Singh, H. and Glass, C. K. (2010).
#' Simple combinations of lineage-determining transcription factors
#' prime cis-regulatory elements required for macrophage and B cell
#' identities. \emph{Molecular Cell} \bold{38}, 576-589.
#' 
#' @export
#' 
run_homer <-
  function (de, k, positions, genome = "hg19",
            subset = function (postmean, lpval, lfsr, rank, quantile)
              lfsr < 0.05,
            homer.exec = "findMotifsGenome.pl",
            out.dir = tempdir(),
            homer.options = "-len 8,10,12 -size 200 -mis 2 -S 25 -p 1 -h",
            verbose = TRUE) {

  # Get the positions if they are not provided.
  if (missing(positions)) {
    feature_names <- rownames(de$postmean)
    out           <- strsplit(feature_names,"_")
    positions     <- data.frame(chr   = sapply(out,"[[",1),
                                start = sapply(out,"[[",2),
                                end   = sapply(out,"[[",3),
                                name  = feature_names,
                                stringsAsFactors = FALSE)
  }

  # Select the differentially expressed positions.
  rows <- select_de_genes(de,k,subset)
  if (verbose)
    cat(sprintf("%d out of %d positions selected\n",
                length(rows),nrow(de$postmean)))

  # Write the selected positions to a BED file.
  pos.file <- file.path(out.dir,"positions.bed")
  if (verbose)
    cat("Writing selected positions to",pos.file,"\n")
  write.table(positions[rows,],pos.file,sep = "\t",quote = FALSE,
              row.names = FALSE,col.names = FALSE)

  # Run the HOMER motif enrichment analysis.
  homer.dir <- file.path(out.dir,"homer")
  homer.command <- paste(homer.exec,pos.file,genome,homer.dir,homer.options)
  if (verbose) {
    cat("Performing HOMER motif enrichment analysis:\n")
    cat(homer.command,"\n")
  }
  system.out <- system(homer.command,ignore.stderr = TRUE,
                       ignore.stdout = TRUE,intern = TRUE)
  res <- read.table(file.path(homer.dir,"knownResults.txt"),
                    sep = "\t",comment.char = "",header = TRUE,
                    check.names = FALSE,stringsAsFactors = FALSE)
  return(res)
}

