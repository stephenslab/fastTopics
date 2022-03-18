# Perform differential accessbility analysis for ATAC-seq regions (peaks),
# and perform TF motif enrichment analysis using HOMER.
select_small_pvals <- function (postmean,lpval,lfsr,rank,quantile)
  lpval > 1
de <- readRDS("DA_regions_topics_noshrinkage_10000iters.rds")

# Select regions with p-value less than 0.1.
regions <- select_de_genes(de,k = 4,subset = select_small_pvals)

# For each topic, perform TF motif enrichment analysis using HOMER
# hypergeometric test.
res <- run_homer(de,k = 4,subset = select_small_pvals)
