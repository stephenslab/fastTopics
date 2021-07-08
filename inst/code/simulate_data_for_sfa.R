library(R.matlab)
set.seed(1)
n   <- 400
m   <- 1000
k   <- 6
dat <- simulate_multinom_gene_data(n,m,k)
X   <- dat$X
L   <- dat$L
writeMat("sim_multinom.mat",X=X,L=L)
