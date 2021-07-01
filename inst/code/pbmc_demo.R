library(Matrix)
# library(fastTopics)
set.seed(1)
data(pbmc_facs)
genes <- pbmc_facs$genes
X     <- pbmc_facs$counts
fit   <- pbmc_facs$fit

# For testing only.
m     <- ncol(X)
i     <- sample(m,2000)
# i   <- c(5018,13171,13978,15685,sample(m,20))
X     <- X[,i]
genes <- genes[i,]
fit$F <- fit$F[i,]

print(system.time(
  out1 <- de_analysis(fit,X,s = rowSums(X) + 1,ns = 1000,nc = 4)))
print(system.time(
  out2 <- de_analysis(fit,X,s = rowSums(X) + 1,ns = 1000,nc = 1)))
k <- 4
dat <- as.data.frame(cbind(low  = out1$low[,k],
                           est  = out1$est[,k],
                           high = out1$high[,k],
                           z    = out1$z[,k]))
rownames(dat) <- genes$symbol
dat <- dat[order(dat$z,decreasing = TRUE),]
# print(dat,digits = 2)
