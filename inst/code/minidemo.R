set.seed(1)
n <- 200
m <- 400
k <- 8
A <- matrix(runif(n*k),n,k)
B <- matrix(runif(m*k),k,m)
X <- A %*% B
