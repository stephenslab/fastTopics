# Perform differential accessbility analysis for ATAC-seq regions (peaks),
# and perform TF motif enrichment analysis using HOMER.
de <- readRDS("DA_regions_topics_noshrinkage_10000iters.rds")

# For each topic, perform TF motif enrichment analysis using HOMER
# hypergeometric test.
select_small_pvals <- function (postmean,lpval,lfsr,rank,quantile)
  lpval > 1
res <- run_homer(de,k = 4,subset = select_small_pvals)
