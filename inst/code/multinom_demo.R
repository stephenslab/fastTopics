# Short script to verify implementation of the differential expression
# (DE) analysis methods applied to data simulated from a multinomial
# topic model.
library(Matrix)
library(ggplot2)
library(cowplot)

# Simulate data.
set.seed(1)
n   <- 400
m   <- 1000
k   <- 4
dat <- simulate_multinom_gene_data(n,m,k,sparse = TRUE)
X   <- dat$X
L   <- dat$L

# Fit a Poisson model (approximating a binomial model) to each gene
# (row) j, and compute the log-fold change statistics.
fit <- init_poisson_nmf(X,L = L,init.method = "random")
de1 <- de_analysis(fit,X,fit.method = "glm")
de2 <- de_analysis(fit,X,fit.method = "scd")
de3 <- de_analysis(fit,X,fit.method = "em")

# Compare the glm and scd estimates of the model parameters.
plot(de1$F + 1e-4,de2$F + 1e-4,pch = 4,cex = 0.5,log = "xy",xlab = "glm",
     ylab = "scd")
abline(a = 0,b = 1,col = "magenta",lty = "dotted")

# Compare the glm and EM estimates of the model parameters.
plot(de1$F + 1e-4,de3$F + 1e-4,pch = 4,cex = 0.5,log = "xy",xlab = "glm",
     ylab = "em")
abline(a = 0,b = 1,col = "magenta",lty = "dotted")

# Compare the scd estimates against the probabilities used to simulate
# the data.
plot(dat$F + 1e-4,de2$F + 1e-4,pch = 4,cex = 0.5,log = "xy",xlab = "true",
     ylab = "estimated")
abline(a = 0,b = 1,col = "magenta",lty = "dotted")

# Here we show that the z-score varies (as expected) with the log-fold
# change estimate and the average expression level.
pdat <- data.frame(f0       = rep(de3$f0,4),
                   postmean = as.vector(de3$postmean),
                   z        = as.vector(de3$z))
print(ggplot(pdat,aes(x = f0,y = postmean,fill = z)) +
  geom_point(size = 2,shape = 21,color = "white") +
  geom_abline(intercept = 0,slope = 0,color = "black",linetype = "dotted") +
  scale_x_continuous(trans = "log10") +
  scale_fill_gradient2(low = "darkblue",mid = "lightskyblue",
                       high = "orangered",midpoint = 0) +
  labs(x = "average expression",y = "log-fold change",fill = "z-score") +
  theme_cowplot(12))
