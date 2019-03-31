# Small script to test that fitpoismix works on toy data sets: one
# encoded as a dense matrix, the other encoded as a sparse matrix.
library(Matrix)
library(quadprog)
library(Rcpp)
source("../code/misc.R")
source("../code/activeset.R")
source("../code/poismix.R")
set.seed(1)

# Create the data set.
U <- matrix(runif(18),6,3) > 0.5
L <- matrix(round(100*runif(18)),6,3)
L <- L * U
w <- c(1,5,100,1,2,0)
x <- c(1,1,1)

# Fit the model.
out1 <- fitpoismix(L,w,c(1,1,1),numiter = 40,qp.solver = "activeset")
out2 <- fitpoismix(L,w,c(1,1,1),numiter = 40,qp.solver = "quadprog")
print(out1$x - out2$x)
print(out1$value - out2$value)

# In this second example, we add a column to L with small entries.
L <- cbind(L,c(1,1,1,0,0,1))
out1 <- fitpoismix(L,w,c(1,0,0,0),numiter = 40,qp.solver = "activeset")
out2 <- fitpoismix(L,w,c(1,1,1,1),numiter = 40,qp.solver = "quadprog")
print(out1$x - out2$x)
print(out1$value - out2$value)
