library(Matrix)
library(fastTopics)
set.seed(1)
data(pbmc_facs)
genes <- pbmc_facs$genes
X     <- pbmc_facs$counts
fit   <- pbmc_facs$fit
set.seed(1)
timing <- system.time(out1 <- de_analysis(fit,X,nc = 1))
print(timing)
set.seed(1)
timing <- system.time(out2 <- de_analysis(fit,X,nc = 4))
print(timing)
testthat::expect_equal(out1,out2,scale = 1,tolerance = 1e-15)
out <- out1
k <- 4
dat <- as.data.frame(cbind(low  = out$low[,k],
                           est  = out$est[,k],
                           high = out$high[,k],
                           z    = out$z[,k]))
rownames(dat) <- paste(genes$symbol,genes$ensembl,sep="_")
dat <- dat[order(dat$z,decreasing = TRUE),]
print(head(dat,n = 4),digits = 4)
print(tail(dat,n = 4),digits = 4)
