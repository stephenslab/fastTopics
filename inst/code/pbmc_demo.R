set.seed(1)
data(pbmc_facs)
X   <- pbmc_facs$counts
fit <- pbmc_facs$fit
out <- de_analysis(fit,X,ns = 1000,nc = 4)
