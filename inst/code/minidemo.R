library(Matrix)
library(NNLM)

# Generate a 300 x 400 data matrix to factorize. Less than 10% of the
# matrix elements should be nonzero.
set.seed(1)
n <- 300
m <- 400
k <- 3
F <- matrix(runif(m*k)/3,m,k)
L <- matrix(runif(n*k)/3,n,k)
X <- matrix(rpois(n*m,L %*% t(F)),n,m)
X <- as(X,"sparseMatrix")
nnzero(X)/(n*m)

# Generate random initial estimates of the factors and loadings.
fit0 <- list(F = matrix(runif(m*k),m,k),
             L = matrix(runif(n*k),n,k))

# Run 60 iterations of the sequential coordinate-wise descent
# algorithm implemented in the NNLM package. Note that nnmf does not
# accept a sparse matrix as input, so we need to provide it witha dense
# matrix instead.
fit1 <- suppressWarnings(
  nnmf(as.matrix(X),k,init = list(W = fit0$L,H = t(fit0$F)),
       method = "scd",loss = "mkl",max.iter = 60,rel.tol = 0, 
       inner.max.iter = 4,trace = 1,verbose = 0))

# Run 60 coordinate-wise updates of the SQP method implemented in the
# fastTopics package. This uses the extrapolation scheme of Ang &
# Gillis (2019) to accelerate the updates.
fit2 <- altsqp(X,fit0,numiter = 60,control = list(extrapolate = 10),
               verbose = FALSE)

# Compare the Poisson log-likelihood at the two solutions; the
# likelihood should be higher at the the altsqp solution.
fit1$F <- t(fit1$H)
fit1$L <- fit1$W
print(loglik.poisson(X,fit1),digits = 14)
print(loglik.poisson(X,fit2),digits = 14)

# Compare the multinomial log-likelihood at the two solutions; again,
# the likelihood should be higher at the altsqp solution.
print(loglik.multinom(X,poisson2multinom(fit1)),digits = 14)
print(loglik.multinom(X,poisson2multinom(fit2)),digits = 14)

# Plot the improvement in the solution over time; the altsqp iterates
# (the solid, orange line) gets much closer to the best solution.
fbest    <- 31041.93896745
fit1$mkl <- n*m*fit1$mkl + sum(X - X*log(X + 1e-16))
plot(fit2$progress$iter,fit2$progress$objective - fbest,
     log = "y",col = "darkorange",type = "l",lwd = 2,xlab = "iteration",
     ylab = "distance from solution")
lines(1:60,fit1$mkl - fbest,col = "darkblue",lwd = 2,lty = "dashed")
