# Verify the implementation of poisson2binom with a small, simulated
# data set.
library(Matrix)
        
# Simulate a 800 x 800 sparse binary matrix from a binomial topic
# model with k = 3 topics.
set.seed(1)
n <- 200
m <- 180

# TO DO: Try with n different from m.
L <- rbind(cbind(rep(1,n),rep(0,n),rep(0,n)),
           cbind(rep(0,n),rep(1,n),rep(0,n)),
           cbind(rep(0,n),rep(0,n),rep(1,n)),
           cbind(runif(n),runif(n),runif(n)))
L <- normalize.rows(L)
F <- cbind(c(rep(0.08,m),rep(0.02,m),rep(0.05,m),rep(0.01,m)),
           c(rep(0.01,m),rep(0.01,m),rep(0.05,m),rep(0.008,m)),
           c(rep(0.08,m),rep(0.00,m),rep(0.05,m),rep(0.1,m)))
P <- L %*% t(F)
n <- nrow(P)
m <- ncol(P)
X <- matrix(rbinom(n*m,1,P),n,m)
print(unique(as.vector(X)))
# TO DO: Try also the case when X is sparse.
# X <- as(X,"dgCMatrix")
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

# ADD COMMENTS HERE.
print(cor(sim$F,fit_binom_em$F))
print(cor(sim$L,fit_binom_em$L))

# ADD COMMENTS HERE.
par(mfrow = c(1,2))
plot(fit_binom_em$F,fit_binom$F,pch = 20)
abline(a = 0,b = 1,col = "magenta",lty = "dashed")
plot(fit_binom_em$L,fit_binom$L,pch = 20)
abline(a = 0,b = 1,col = "magenta",lty = "dashed")

# ADD COMMENTS HERE.
fit_multinom <- poisson2multinom(fit_pois)
par(mfrow = c(1,2))
plot(fit_binom_em$F,fit_multinom$F,pch = 20)
abline(a = 0,b = 1,col = "magenta",lty = "dashed")
plot(fit_binom_em$L,fit_multinom$L,pch = 20)
abline(a = 0,b = 1,col = "magenta",lty = "dashed")
