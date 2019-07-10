library(Matrix)
library(NNLM)

set.seed(1)
n <- 300
m <- 400
k <- 3
A <- matrix(runif(n*k)/3,n,k)
B <- matrix(runif(m*k)/3,k,m)
X <- matrix(rpois(n*m,A %*% B),n,m)
X <- as(X,"sparseMatrix")

F0 <- matrix(runif(m*k),m,k)
L0 <- matrix(runif(n*k),n,k)

fit1 <- nnmf(as.matrix(X),k,init = list(W = L0,H = t(F0)),method = "scd",
             loss = "mkl",
             max.iter = 500,rel.tol = 0,n.threads = 0,inner.max.iter = 4,
             trace = 1,verbose = 2)
fit2 <- altsqp(X,F0,L0,numiter = 100,control = list(extrapolate = 10))
print(cost(X,fit1$W %*% fit1$H,1e-15),digits = 12)
print(fit2$value,digits = 12)
