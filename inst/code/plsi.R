# Small script to verify that the Poisson NMF multiplicative updates
# are equivalent to the pLSI EM updates.

# Simulate a 100 x 200 counts matrix.
set.seed(1)
n   <- 100
m   <- 200
k   <- 3
out <- simulate_count_data(n,m,k)
X   <- out$X
A   <- out$L
B   <- t(out$F)

# Apply the EM updates for pLSI and Poisson NMF in parallel.
N <- rowSums(X)
for (iter in 1:20) {
  out <- get_multinom_from_pnmf(t(B),A)
  L   <- out$L
  F   <- out$F
    
  # Apply the multiplicative (EM) update for L.
  A <- scale.cols(A * tcrossprod(X / (A %*% B),B),1/rowSums(B))

  # Apply the pLSI EM update for L.
  P <- matrix(0,m,k)
  for (i in 1:n) {
    for (j in 1:m)
      P[j,] <- F[j,]*L[i,]/sum(F[j,]*L[i,])
    L[i,] <- (X[i,] %*% P)/N[i]
  }

  # Compare the updated L matrices.
  out <- get_multinom_from_pnmf(t(B),A)
  cat(sprintf("%0.1e ",max(abs(out$L - L))))
  
  # Apply the multiplicative (EM) update for F.
  B <- B * crossprod(A,X / (A %*% B)) / colSums(A)
  
  # Apply the pLSI EM update for F.
  P <- matrix(0,n,k)
  for (j in 1:m) {
    for (i in 1:n)
      P[i,] <- F[j,]*L[i,]/sum(F[j,]*L[i,])
    F[j,] <- X[,j] %*% P
  }
  F <- normalize.cols(F)

  # Compare the updated F matrices
  out <- get_multinom_from_pnmf(t(B),A)
  cat(sprintf("%0.1e\n",max(abs(out$F - F))))
}
