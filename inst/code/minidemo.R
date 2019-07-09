library(NNLM)

set.seed(1)
n <- 100
m <- 200
k <- 3
A <- matrix(runif(n*k),n,k)
B <- matrix(runif(m*k),k,m)
X <- A %*% B

F0 <- matrix(runif(m*k),m,k)
L0 <- matrix(runif(n*k),n,k)

fit1 <- betanmf(X,L0,t(F0))
fit2 <- nnmf(X,k,init = list(W = L0,H = t(F0)),method = "scd",loss = "mkl",
             max.iter = 1000,rel.tol = 0,n.threads = 0,inner.max.iter = 4,
             trace = 1,verbose = 2)
fit3 <- altsqp(X,F0,L0)

print(cost(X,fit1$A %*% fit1$B,1e-15),digits = 12)
print(cost(X,fit2$W %*% fit2$H,1e-15),digits = 12)
print(fit3$value,digits = 12)
