#' @title Run HOMER for Motif Enrichment
#' 
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
#' @export
#' 
run_homer <- function(pos.file,
                      genome,
                      homer.exec = "findMotifsGenome.pl",
                      out.dir = "out",
                      region.size = 200,
                      motif.length = c(8,10,12),
                      optimize.count = 2,
                      homer.options = "-h -S 25 -p 1") {

  # if(!dir.exists(out.dir))
  #   dir.create(out.dir,recursive = TRUE)
  #   homer.path <- normalizePath(homer.path)
  #   if(!file_test('-x', homer.path)){
  #     stop(homer.path, " does not exist or is not executable!")
  #   }

  system(paste(homer.exec,pos.file,genome,out.dir,
    '-len', paste0(motif.length, collapse = ','),
    '-size', region.size,
    '-mis', optimize.count,
    '-S', n.motifs,
  ))
  res <-
    read.csv(file.path(out.dir, "knownResults.txt"), sep="\t", header=TRUE, check.names=F, stringsAsFactors=F)
    return(res)
}

