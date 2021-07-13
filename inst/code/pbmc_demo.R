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
out1 <- de_analysis(fit,X,shrink.method = "none",control = list(nc = 4))
set.seed(1)
out2 <- de_analysis(fit,X,shrink.method = "ash",control = list(nc = 4))

# Plot the distribution of MCMC acceptance rates.
hist(out2$ar,n = 64)

# Compare the LFC estimates with and without shrinkage.
k <- 4
pdat1 <- data.frame(b1 = out1$est[,k],b2 = out2$est[,k])
p1 <- ggplot(pdat1,aes(x = b1,y = b2)) +
  geom_point(shape = 21,color = "white",fill = "black",na.rm = TRUE) +
  geom_abline(intercept = 0,slope = 1,color = "magenta",linetype = "dotted") +
  labs(x = "MLE estimate",y = "stabilized MLE estimate") +
  theme_cowplot(font_size = 10)

# Create a volcano plot to visualize the DE results for topic k = 4.
pdat2 <- data.frame(lfc  = out2$est[,k],
                    z    = abs(out2$z[,k]),
                    f0   = log10(out2$f0),
                    gene = genes$symbol,
                    stringsAsFactors = FALSE)
rows <- which(!with(pdat2,lfc > 1 & z > 8))
pdat2[rows,"gene"] <- ""
p2 <- ggplot(pdat2,aes(x = lfc,y = z,fill = f0,label = gene)) +
  geom_point(color = "white",stroke = 0.3,shape = 21,
             na.rm = TRUE) +
  geom_text_repel(color = "darkgray",size = 2.25,fontface = "italic",
                  segment.color = "darkgray",segment.size = 0.25,
                  min.segment.length = 0,max.overlaps = Inf,
                  na.rm = TRUE) +
  scale_y_continuous(trans = "sqrt",limits = c(0,80),
                     breaks = c(0,1,2,5,10,20,50,100)) +
  scale_fill_gradient2(low = "deepskyblue",mid = "gold",high = "orangered",
                       na.value = "gainsboro",
                       midpoint = mean(range(pdat2$f0,na.rm = TRUE))) +
  labs(x = "log-fold change",y = "|z-score|",fill = "log10(f0)") +
  theme_cowplot(font_size = 10)
