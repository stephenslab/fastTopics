# Short script to verify implementation of the differential count
# analysis methods applied to data simulated from a multinomial topic
# model.
library(Matrix)
library(ggplot2)
library(cowplot)

# Simulate data.
set.seed(1)
n   <- 800
m   <- 1000
k   <- 4
dat <- simulate_multinom_gene_data(n,m,k,sparse = TRUE)
X   <- dat$X
L   <- dat$L

# Fit a Poisson model (approximating a binomial model) to each
# combination of gene j and topic k, and compute the log-fold change
# statistics.
fit <- init_poisson_nmf(X,L = L,init.method = "random")
out <- diff_count_analysis(fit,X)

# For the selected topic, compare the f1 estimates against the
# probabilities used to simulate the data.
i <- 1
plot(dat$F[,i],out$F1[,i],pch = 4,cex = 0.5,log = "xy",
     xlab = "true f1",ylab = "estimated f1")
abline(a = 0,b = 1,col = "magenta",lty = "dotted")

# For the selected topic, compare f0 estimates against the probabities
# used to simulate the data. Since there is no "f0" used to simulate
# the data, here we approximate f0 by taking the average of the
# Poisson rates across all topics other than topic i.
plot(rowMeans(dat$F[,-i]),out$F0[,i],pch = 4,cex = 0.5,log = "xy",
     xlab = "true f0 (approx)",ylab = "estimated f0")
abline(a = 0,b = 1,col = "magenta",lty = "dotted")

# Here we show that the Z-score varies (predictably) with the log-fold
# change estimate and the average expression level.
pdat <- data.frame(x = out$colmeans,beta = out$beta[,i],z = out$Z[,i])
print(ggplot(pdat,aes(x = x,y = beta,fill = z)) +
  geom_point(size = 2,shape = 21,color = "white") +
  geom_abline(intercept = 0,slope = 0,color = "black",linetype = "dotted") +
  scale_x_continuous(trans = "log10") +
  scale_fill_gradient2(low = "darkblue",mid = "lightskyblue",
                       high = "orangered",midpoint = 0) +
  labs(x = "average expression",y = "log-fold change",fill = "z-score") +
  theme_cowplot(12))

# Show a traditional volcano plot, in which log-fold change is shown
# on the x-axis, and the z-score is shown on the y-axis. To illustrate
# the impact of overall gene expression level on the z-scores, the
# (log) average expression level is shown by a colour gradient.
print(ggplot(pdat,aes(x = beta,y = abs(z),fill = log10(x))) +
  geom_point(size = 1.5,shape = 21,color = "white") +
  labs(x = "log-fold change",y = "|z-score|",fill = "log10(avg.exp.)") +
  scale_fill_gradient2(low = "skyblue",mid = "gold",
                       high = "orangered",midpoint = 0) +
  theme_cowplot(12))
