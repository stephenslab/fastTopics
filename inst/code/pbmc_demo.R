data(pbmc_facs)
X   <- pbmc_facs$counts
fit <- pbmc_facs$fit
out <- de_analysis(fit,X,nc = 4)
