library(Matrix)
# library(fastTopics)
set.seed(1)
data(pbmc_facs)
genes <- pbmc_facs$genes
X     <- pbmc_facs$counts
fit   <- pbmc_facs$fit

# For testing only.
# m     <- ncol(X)
# i     <- sample(m,2000)
# i     <- c(5018,13171,13978,15685,sample(m,20))
# X     <- X[,i]
# genes <- genes[i,]
# fit$F <- fit$F[i,]

print(system.time(
  out <- de_analysis(fit,X,s = rowSums(X) + 1,ns = 100,nc = 4)))
k <- 4
dat <- as.data.frame(cbind(low  = out$low[,k],
                           est  = out$est[,k],
                           high = out$high[,k],
                           z    = out$z[,k]))
rownames(dat) <- paste(genes$symbol,genes$ensembl,sep="_")
dat <- dat[order(dat$z,decreasing = TRUE),]
print(head(dat,n = 10),digits = 2)
print(tail(dat,n = 10),digits = 2)
