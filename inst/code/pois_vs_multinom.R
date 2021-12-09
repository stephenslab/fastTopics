# A short script to verify the fit_poisson_models computations against
# fit_topic_model.
library(Matrix)

# Simulate data.
set.seed(1)
n   <- 120
m   <- 1000
k   <- 4
dat <- simulate_multinom_gene_data(n,m,k,sparse = TRUE)
X   <- dat$X
L   <- dat$L

# Fit a multinomial topic model, with k = 4 topics.
fit <- fit_topic_model(X,k = 4)

# Ensure that none of the topic proportions are exactly zero or
# exactly one.
L <- fit$L
L <- pmax(L,1e-8)
L <- pmin(L,1 - 1e-8)

s   <- rowSums(X)
out <- add_pseudocounts(X,s*L,0.01)
X   <- out$X
L   <- out$L

# For each column j of the counts matrix, compute MLEs of the
# parameters in the Poisson glm, x ~ Poisson(u), in which the
# Poisson rates are u = sum(L*f), and f = F[j,].
F <- fit_poisson_models(X,L,"scd",1e-8,20,1e-8,1)
F <- pmax(F,1e-8)

# Compare the estimates obtained by computing MLEs under the
# multinomial topic model against the estimates obtained by running
# fit_poisson_models.
plot(fit$F + 1e-6,F + 1e-6,pch = 4,cex = 0.5,log = "xy")
abline(a = 0,b = 1,col = "skyblue",lty = "dotted")
