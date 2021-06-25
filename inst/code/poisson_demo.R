# Short script to verify implementation of the differential expression
# (DE) analysis methods applied to data simulated from a Poisson NMF
# model.
library(Matrix)
library(ggplot2)
library(cowplot)

# Simulate data.
set.seed(1)
n   <- 800
m   <- 1000
k   <- 4
s   <- 10^runif(n,-1,1)
dat <- simulate_poisson_gene_data(n,m,k,s)
X   <- dat$X
L   <- dat$L
Y   <- as(X,"dgCMatrix")
mu  <- colMeans(X)
f0  <- colSums(X)/sum(s)

# Add "pseudocounts" to the data.
out <- add_pseudocounts(X,s*L,0.01)
X <- out$X
L <- out$L

# Fit a Poisson model for each gene.
F1 <- fit_poisson_models(X,L,method = "glm")
F2 <- fit_poisson_models(X,L,method = "scd",nc = 4)
print(range(F1 - F2))

# Compare the estimates against the Poisson rates used to simulate the
# data.
e <- 1e-4
i <- 1
plot(dat$F + e,F1 + e,pch = 20,log = "xy",xlab = "true frequency",
     ylab = "estimated frequency")
abline(a = 0,b = 1,col = "dodgerblue",lty = "dotted")

# Compute the log-fold change statistics for each gene j and topic k.
out <- compute_lfc_stats(X,F1,L,f0)

# Here we show that the z-score varies, as expected, with the log-fold
# change estimate and the average expression level.
pdat <- data.frame(x = mu,lfc = out$est[,1],z = out$z[,1])
print(ggplot(pdat,aes(x = x,y = lfc,fill = z)) +
  geom_point(size = 2,shape = 21,color = "white") +
  geom_abline(intercept = 0,slope = 0,color = "black",linetype = "dotted") +
  scale_x_continuous(trans = "log10") +
  scale_fill_gradient2(low = "darkblue",mid = "skyblue",
                       high = "orangered",midpoint = 0) +
  labs(x = "average expression",y = "log-fold change",fill = "z-score") +
  theme_cowplot(font_size = 12))

# Create a volcano plot in which log-fold change is shown on the
# x-axis and the z-score is shown on the y-axis. To illustrate the
# impact of (mean) gene expression level on the z-scores, the (log)
# average expression level is shown by a colour gradient.
print(ggplot(pdat,aes(x = lfc,y = abs(z),fill = log10(x))) +
  geom_point(size = 2,shape = 21,color = "white") +
  labs(x = "log-fold change",y = "|z-score|",fill = "log10(mean)") +
  scale_y_continuous(trans = "sqrt") +
  scale_fill_gradient2(low = "skyblue",mid = "gold",high = "orangered",
                       midpoint = 0) +
  theme_cowplot(font_size = 12))
