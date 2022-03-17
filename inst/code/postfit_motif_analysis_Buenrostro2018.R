# Perform differential accessbility analysis for ATAC-seq regions (peaks),
# and perform TF motif enrichment analysis using HOMER.
de <- readRDS("DA_regions_topics_noshrinkage_10000iters.rds")

# Select regions with p-value less than 0.01.
regions <-
  select_de_genes(de,k = 4,
                  subset = function (postmean,lpval,lfsr,rank,quantile)
                             lpval > 1)

# For each topic, perform TF motif enrichment analysis using HOMER
# hypergeometric test.
homer_res        <- vector("list", ncol(DA_res$z))
names(homer_res) <- colnames(DA_res$z)
for (k in 1:ncol(DA_res$z)) {
  homer_res[[k]] <-
    run_homer(selected_regions$filenames[k],
              genome = "hg19",
              homer.path = homerpath,
              use.hypergeometric = TRUE,
              out.dir = "/Users/pcarbo/homer/findMotifsGenome.pl",
              n.cores = 2)
}
