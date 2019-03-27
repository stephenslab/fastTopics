# Small script to test that fitpoismix works on toy data sets: one
# encoded as a dense matrix, the other encoded as a sparse matrix.
library(Matrix)
library(quadprog)
source("../code/misc.R")
source("../code/altsqp.R")
set.seed(1)

# Example 1: dense matrix.
U <- matrix(runif(18),6,3) > 0.5
L <- matrix(round(100*runif(18)),6,3)
L <- L * U
w <- c(1,5,100,1,2,0)
x <- c(1,1,1)
out1 <- fitpoismix(L,w,x,numiter = 20)

# Example 2: sparse matrix.
L    <- Matrix(L,sparse = TRUE)
out2 <- fitpoismix(L,w,x,numiter = 20)

