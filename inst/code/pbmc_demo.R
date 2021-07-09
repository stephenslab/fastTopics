library(Matrix)
set.seed(1)
data(pbmc_facs)
genes <- pbmc_facs$genes
X     <- pbmc_facs$counts
fit   <- pbmc_facs$fit
# set.seed(1)
# timing <- system.time(out1 <- de_analysis(fit,X,control = list(nc = 1)))
# print(timing)
# set.seed(1)
# timing <- system.time(out2 <- de_analysis(fit,X,control = list(nc = 4)))
# print(timing)
# testthat::expect_equal(out1,out2,scale = 1,tolerance = 1e-15)
# out <- out2
out <- de_analysis(fit,X,lfc.stat = "de",control = list(ns = 1000,nc = 4))
k <- 4
dat <- as.data.frame(cbind(lower = out$lower[,k],
                           est   = out$est[,k],
                           mean  = out$est[,k],
                           upper = out$upper[,k],
                           z     = out$z[,k],
                           lf0   = log10(out$f0)))
rownames(dat) <- paste(genes$symbol,genes$ensembl,sep="_")
dat <- dat[order(dat$z,decreasing = TRUE),]
print(head(dat,n = 8),digits = 4)
print(tail(dat,n = 4),digits = 4)
hist(out$ar,n = 64)
