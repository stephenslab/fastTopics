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
set.seed(1)
out1 <- de_analysis(fit,X,shrink.method = "none",control = list(nc = 4))
set.seed(1)
out2 <- de_analysis(fit,X,control = list(nc = 4))
out <- out2
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
                   gene = genes$symbol,
                   stringsAsFactors = FALSE)
rows <- which(!with(pdat,lfc > 1 & z > 8))
pdat[rows,"gene"] <- ""
p1 <- ggplot(pdat,aes(x = lfc,y = z,fill = f0,label = gene)) +
  geom_point(color = "white",stroke = 0.3,shape = 21,
             na.rm = TRUE) +
  geom_text_repel(color = "darkgray",size = 2.25,fontface = "italic",
                  segment.color = "darkgray",segment.size = 0.25,
                  min.segment.length = 0,max.overlaps = Inf,
                  na.rm = TRUE) +
  scale_y_continuous(trans = "sqrt",
                     breaks = c(0,1,2,5,10,20,50,100)) +
  scale_fill_gradient2(low = "deepskyblue",mid = "gold",high = "orangered",
                       na.value = "gainsboro",
                       midpoint = mean(range(pdat$f0,na.rm = TRUE))) +
  labs(x = "log-fold change",y = "|z-score|",fill = "log10(f0)") +
  theme_cowplot(font_size = 10)
pdat <- data.frame(b1 = out1$est[,k],b2 = out2$est[,k])
p2 <- ggplot(pdat,aes(x = b1,y = b2)) +
  geom_point(shape = 21,color = "white",fill = "black",na.rm = TRUE) +
  geom_abline(intercept = 0,slope = 1,color = "magenta",linetype = "dotted") +
  labs(x = "original estimate",y = "stabilized estimate") +
  theme_cowplot(font_size = 10)
