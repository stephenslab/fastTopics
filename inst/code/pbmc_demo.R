library(ggplot2)
library(ggrepel)
library(cowplot)
library(Matrix)

# Load the data.
set.seed(1)
data(pbmc_facs)
genes <- pbmc_facs$genes
X     <- pbmc_facs$counts
fit   <- pbmc_facs$fit

# Perform the differential expression analysis, with and without
# shrinkage.
set.seed(1)
out1 <- de_analysis(fit,X,shrink.method = "none",
                    control = list(nc = 4,nsplit = 400))
set.seed(1)
out2 <- de_analysis(fit,X,shrink.method = "ash",
                    control = list(nc = 4,nsplit = 400))

# Plot the distribution of MCMC acceptance rates.
hist(out2$ar,n = 64)

# Compare the LFC estimates with and without shrinkage.
pdat <- data.frame(postmean1 = as.vector(out1$postmean),
                   postmean2 = as.vector(out2$postmean))
p1 <- ggplot(pdat,aes(x = postmean1,y = postmean2)) +
  geom_point(shape = 21,color = "white",fill = "darkblue",na.rm = TRUE) +
  geom_abline(intercept = 0,slope = 1,color = "magenta",linetype = "dotted") +
  labs(x = "posterior mean estimate",y = "stabilized posterior estimate") +
  theme_cowplot(font_size = 10)

# Create a volcano plot to visualize the DE results for topic k = 4.
k <- "k4"
pdat <- data.frame(gene     = genes$symbol,
                   postmean = out2$postmean[,k],
                   z        = pmin(40,abs(out2$z[,k])),
                   lfsr     = cut(out2$lfsr[,k],c(-1,0.001,0.01,0.05,Inf)),
                   stringsAsFactors = FALSE)
rows <- which(with(pdat,!(postmean > 3 | (postmean > 0 & z > 10))))
pdat[rows,"gene"] <- ""
p2 <- ggplot(pdat,aes(x = postmean,y = z,fill = lfsr,label = gene)) +
  geom_point(color = "white",stroke = 0.3,shape = 21,
             na.rm = TRUE) +
  geom_text_repel(color = "darkgray",size = 2.25,fontface = "italic",
                  segment.color = "darkgray",segment.size = 0.25,
                  min.segment.length = 0,max.overlaps = Inf,
                  na.rm = TRUE) +
  scale_y_continuous(trans = "sqrt",breaks = c(0,1,2,5,10,20,50)) +
  scale_fill_manual(values = c("deepskyblue","gold","orange","tomato"),
                    na.value = "gainsboro") +
  labs(x = "log-fold change",y = "|z-score|") +
  theme_cowplot(font_size = 10)
