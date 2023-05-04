# Verify the implementation of poisson2binom with a small, simulated
# data set.
library(Matrix)
        
# Simulate a 800 x 800 sparse binary matrix from a binomial topic
# model with k = 3 topics.
set.seed(1)
n <- 200
m <- 200
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
# TO DO: Try also the case when M is sparse.
# X <- as(X,"dgCMatrix")
sim <- list(L = L,F = F,X = X)

# Fit a Poisson non-negative matrix factorization to the binomial
# data.  To simplify comparison with the "true" factorization---that
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
