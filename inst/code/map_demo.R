# TO DO: Explain here briefly what this script is for.
library(Matrix)
set.seed(1)
n <- 40
m <- 80
k <- 4
X <- simulate_count_data(n,m,k = k,sparse = TRUE)$X
fit0 <- fit_topic_model(X,k = k,numiter.main = 10,numiter.refine = 10)
fit1 <- fit_topic_model_map(X,fit0,numiter = 100,
                            alpha.factors = matrix(1/m,m,k),
                            alpha.loadings = matrix(1/k,n,k))
cat("\n")
print(range(diff(fit1$logposterior)))
y <- max(fit1$logposterior) - fit1$logposterior + 0.01 
plot(1:100,y,type = "l",log = "y",lwd = 2,xlab = "iteration",
     ylab = "distance to best log-posterior")
