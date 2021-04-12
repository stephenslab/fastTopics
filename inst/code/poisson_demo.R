# Short script to verify implementation of the differential count
# analysis methods applied to data simulated from a Poisson NMF model.
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

# Fit a Poisson model for each gene.
F1 <- fit_poisson_models(X,L,s,method = "glm")
F2 <- fit_poisson_models(X,L,s,method = "scd",nc = 4)
print(range(F1 - F2))

# Compare the estimates against the Poisson rates used to simulate the
# data.
e <- 1e-3
i <- 1
plot(dat$F + e,F1 + e,pch = 20,log = "xy",xlab = "true F",ylab = "estimated F")
abline(a = 0,b = 1,col = "dodgerblue",lty = "dotted")

# Compute the log-fold change and Z-score for each gene j and topic k.
out1 <- compute_univar_poisson_zscores(X,L,out.em$F0,out.em$F1,s)
out2 <- compute_univar_poisson_zscores_fast(X,L,out.em$F0,out.em$F1,s)
out3 <- compute_univar_poisson_zscores_fast(Y,L,out.em$F0,out.em$F1,s)
print(max(abs(out1$beta - out2$beta)))
print(max(abs(out1$se   - out2$se)/out2$se,na.rm = TRUE))
print(max(abs(out1$Z    - out2$Z)))
print(max(abs(out1$pval - out2$pval)))
print(max(abs(out1$beta - out3$beta)))
print(max(abs(out1$se   - out3$se)/out1$se,na.rm = TRUE))
print(max(abs(out1$Z    - out3$Z)))
print(max(abs(out1$pval - out3$pval)))

# Here we show that the Z-score varies (predictably) with the log-fold
# change estimate and the average expression level.
pdat <- data.frame(x = colMeans(X),beta = out1$beta[,i],z = out1$Z[,i])
pdat <- subset(pdat,beta > -5)
print(ggplot(pdat,aes(x = x,y = beta,fill = z)) +
  geom_point(size = 2,shape = 21,color = "white") +
  geom_abline(intercept = 0,slope = 0,color = "gray",linetype = "dotted") +
  scale_x_continuous(trans = "log10") +
  scale_fill_gradient2(low = "darkblue",mid = "skyblue",
                       high = "orangered",midpoint = 0) +
  labs(x = "average expression",y = "log-fold change",fill = "z-score") +
  theme_cowplot(12))

# Show a traditional volcano plot, in which log-fold change is shown
# on the x-axis, and the z-score is shown on the y-axis. To illustrate
# the impact of overall gene expression level on the z-scores, the
# (log) average expression level is shown by a colour gradient.
print(ggplot(pdat,aes(x = beta,y = abs(z),fill = log10(x))) +
  geom_point(size = 2,shape = 21,color = "white") +
  labs(x = "log-fold change",y = "|z-score|",fill = "log10(avg.exp.)") +
  scale_fill_gradient2(low = "skyblue",mid = "gold",
                       high = "orangered",midpoint = 0) +
  theme_cowplot(12))
