# TO DO: Explain here what this script does, and how to use it.

# Simulate a 100 x 200 counts matrix.
set.seed(1)
n   <- 100
m   <- 200
k   <- 3
out <- simulate_count_data(n,m,k)
X   <- out$X
F   <- out$F
L   <- out$L

# Add pseudocounts.
a <- 1.1
b <- 1
X <- rbind(X,matrix(a - 1,k,m))
L <- rbind(L,diag(k))
# X <- cbind(X,matrix(b - 1,n+k,k))
# F <- rbind(F,1e-6 * diag(k))

# Fit a multinomial topic model, with k = 3.
fit <- init_poisson_nmf(X,F = F,L = L)
fit <- fit_poisson_nmf(X,fit0 = fit,numiter = 400,
                       update.loadings = 1:n,
                       control = list(extrapolate = TRUE))
fit.multinom <- poisson2multinom(fit)

# Apply the pLSI EM update for L.
X <- X[1:n,1:m]
F <- fit.multinom$F[1:m,]
L <- fit.multinom$L[1:n,]
P <- matrix(0,m,k)
for (i in 1:n) {
  for (j in 1:m)
    P[j,] <- F[j,]*L[i,]/sum(F[j,]*L[i,])
  L[i,] <- X[i,] %*% P + b - 1
}
L <- normalize.rows(L)
print(range(L - fit.multinom$L[1:n,]))

# Apply the pLSI EM update for F.
P <- matrix(0,n,k)
for (j in 1:m) {
  for (i in 1:n)
    P[i,] <- F[j,]*L[i,]/sum(F[j,]*L[i,])
  F[j,] <- X[,j] %*% P + a - 1
}
F <- normalize.cols(F)
print(range(F - fit.multinom$F[1:m,]))
