library(ggplot2)
library(ggrepel)
library(cowplot)
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
pdat <- data.frame(lfc  = out$est[,k],
                   z    = abs(out$z[,k]),
                   f0   = log10(out$f0),
                   gene = pbmc_facs$genes$symbol,
                   stringsAsFactors = FALSE)
rows <- which(!with(pdat,lfc > 1 & z > 4))
pdat[rows,"gene"] <- ""
p <- ggplot(pdat,aes(x = lfc,y = z,fill = f0,label = gene)) +
  geom_point(color = "white",stroke = 0.3,shape = 21,
             na.rm = TRUE) +
  geom_text_repel(color = "darkgray",size = 2.25,fontface = "italic",
                  segment.color = "darkgray",segment.size = 0.25,
                  min.segment.length = 0,max.overlaps = Inf,
                  na.rm = TRUE) +
  scale_y_continuous(trans = "sqrt",limits = c(0,20),
                     breaks = c(0,1,2,5,10,20)) +
  scale_fill_gradient2(low = "deepskyblue",mid = "gold",high = "orangered",
                       na.value = "gainsboro",
                       midpoint = mean(range(pdat$f0))) +
  labs(x = "log-fold change",y = "|z-score|",fill = "log10(f0)") +
  theme_cowplot(font_size = 10)
