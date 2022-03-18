#' @title Perform HOMER Motif Enrichment Analysis Using DE Positions
#'
#' @description http://homer.ucsd.edu/homer/
#' See http://homer.ucsd.edu/homer/ngs/peakMotifs.html for more details.
#' 
#' @param regions.file BED file of the peak regions.
#' @param genome Genome assembly; e.g., "hg19"
#' @param homer.path Directory path to 'findMotifsGenome.pl' excutable file.
#' @param out.dir Output directory.
#' 
#' @param region.size Size of the regions used for motif finding (default=200)
#' @param motif.length Length of motifs to be found (default=8,10,12).
#' HOMER will find motifs of each size separately and then combine the results at the end.
#' The length of time it takes to find motifs increases greatly with increasing size.
#' @param optimize.count Number of mismatches allowed in global optimization phase (default=2)
#' @param n.motifs Number of motifs to find (default=25)
#' @param use.hypergeometric Whether to use hypergeometric distribution to score motif enrichment
#' By default, findMotifsGenome.pl uses the binomial distribution to score motifs, which is faster.
#' This works well when the number of background sequences greatly out number the target sequences.
#' If the number of background sequences is smaller than target sequences, it is recommended to use the hypergeometric distribution.
#' @param background BED file of the background regions
#' 
#' @return A data frame of motif enrichment results loaed from HOMER's output 'knownResults.txt'
#'
#' @importFrom utils read.table
#' @importFrom utils write.table
#' 
#' @export
#' 
run_homer <-
  function (de, k, positions, genome = "hg19",
            subset = function (postmean, lpval, lfsr, rank, quantile)
              lfsr < 0.05,
            homer.exec = "findMotifsGenome.pl",
            out.dir = tempdir(),
            homer.options = "-len 8,10,12 -size 200 -mis 2 -h -S 25 -p 1",
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
  write.table(positions[rows,],pos.file,sep = " ",quote = FALSE,
              row.names = FALSE,col.names = FALSE)

  homer.command <- paste(homer.exec,pos.file,genome,out.dir,homer.options)
  if (verbose) {
    cat("Performing HOMER motif enrichment analysis:\n")
    cat(homer.command,"\n")
  }
  out <- system(homer.command)
  res <-
    read.table(file.path(out.dir, "knownResults.txt"),sep="\t", header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
    return(res)
}

