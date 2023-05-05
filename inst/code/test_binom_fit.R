# Verify the implementation of poisson2binom with a small, simulated
# data set.
library(Matrix)
library(ggplot2)
library(cowplot)
        
# Simulate a 800 x 800 sparse binary matrix from a binomial topic
# model with k = 3 topics.
set.seed(1)
n <- 1000
m <- 100
L <- rbind(cbind(rep(1,n),rep(0,n),rep(0,n)),
           cbind(rep(0,n),rep(1,n),rep(0,n)),
           cbind(rep(0,n),rep(0,n),rep(1,n)),
           cbind(runif(n),runif(n),runif(n)))
L <- normalize.rows(L)
F <- rbind(diag(3)/3,
           cbind(c(rep(0.1,m),rep(0.1,m),rep(0.05,m),rep(0.0,m),rep(0.0,m)),
                 c(rep(0.0,m),rep(0.0,m),rep(0.05,m),rep(0.1,m),rep(0.0,m)),
                 c(rep(0.1,m),rep(0.0,m),rep(0.05,m),rep(0.0,m),rep(0.1,m))))
P <- L %*% t(F)
n <- nrow(P)
m <- ncol(P)
X <- matrix(rbinom(n*m,1,P),n,m)
rownames(X) <- paste0("s",1:n)
rownames(L) <- paste0("s",1:n)
colnames(X) <- paste0("m",1:m)
rownames(F) <- paste0("m",1:m)
colnames(L) <- paste0("k",1:3)
colnames(F) <- paste0("k",1:3)
X <- as(X,"dgCMatrix")
sim <- list(L = L,F = F,X = X)

# Fit a Poisson non-negative matrix factorization to the binomial
# data. To simplify comparison with the "true" factorization---that
# is, the L and F used to simulate the data--the factorization is
# initialized to the true parameter values.
fit_pois <- init_poisson_nmf(X,L = L,F = F)
fit_pois <- fit_poisson_nmf(X,fit0 = fit_pois,control = list(extrapolate = TRUE))

# Convert the Poisson NMF to a binomial topic model without any EM
# updates to refine the fit.
fit_binom <- poisson2binom(X,fit_pois,numem = 0)

# Perform the conversion a second time, this time with some EM updates
# to refine the fit.
fit_binom_em <- poisson2binom(X,fit_pois,numem = 20)

# A quick sanity check that the estimates align well with the ground
# truth.
print(cor(sim$F,fit_binom_em$F))
print(cor(sim$L,fit_binom_em$L))

# The Poisson NMF estimates should be close to the binomial topic
# model estimates here because the binomial probabilities are small
# and the sample size is reasonably large.
par(mfrow = c(1,2))
plot(fit_binom_em$F,fit_binom$F,pch = 20)
abline(a = 0,b = 1,col = "magenta",lty = "dashed")
plot(fit_binom_em$L,fit_binom$L,pch = 20)
abline(a = 0,b = 1,col = "magenta",lty = "dashed")

# By contrast, the multinomial topic model estimates of the word
# proportions (F) are very wrong, and, more worryingly, the topic
# proportions (L) are also off.
fit_multinom <- poisson2multinom(fit_pois)
par(mfrow = c(1,2))
plot(fit_binom_em$F,fit_multinom$F,pch = 20)
abline(a = 0,b = 1,col = "magenta",lty = "dashed")
plot(fit_binom_em$L,fit_multinom$L,pch = 20)
abline(a = 0,b = 1,col = "magenta",lty = "dashed")

# Perform the GoM differential expression analysis using the binomial topic model. 
de_binom <- de_analysis(fit_binom_em,X,control = list(nc = 4,ns = 1e4))

# Verify that the F estimates from de_analysis are close to the ground truth.
plot(fit_binom_em$F,de_binom$F,pch = 20)
abline(a = 0,b = 1,col = "magenta",lty = "dashed")

# Plot posterior mean l.e. LFC vs. z-score, coloring the points by the
# "true" change in order to show visually that the GoM DE analysis is
# working.
df <- sim$F - apply(sim$F,1,function (x) sort(x,decreasing = TRUE)[2])
p <- qplot(x = de_binom$postmean,y = de_binom$z,color = df) +
  scale_color_gradient2(low = "darkblue",mid = "deepskyblue",high = "orangered",
                        midpoint = 0) +
  theme_cowplot(font_size = 12)
print(p)

# Perform a second GoM differential expression analysis using the
# multinomial topic model.
de_multinom <- de_analysis(fit_multinom,X,control = list(nc = 4,ns = 1e4))

# Plot posterior mean l.e. LFC vs. z-score, coloring the points by the
# "true" change in order to show visually that the GoM DE analysis is
# working.
p2 <- qplot(x = de_multinom$postmean,y = de_multinom$z,color = df) +
  scale_color_gradient2(low = "darkblue",mid = "deepskyblue",high = "orangered",
                        midpoint = 0) +
  theme_cowplot(font_size = 12)
print(p2)
