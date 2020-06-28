# TO DO: Explain here what this script is for, and how to use it.

# Simulate data.
set.seed(1)
n   <- 800
m   <- 1000
k   <- 4
s   <- 10^runif(n,-1,1)
dat <- simulate_poisson_gene_data(n,m,k,s)
X   <- dat$X
L   <- dat$L

# Fit a Poisson non-negative matrix factorization to the data.
F0 <- matrix(0,m,k)
F1 <- matrix(0,m,k)
for (i in 1:m)
  for (j in 1:k) {
    out     <- fit_poisson_optim(X[,i],s,L[,j])
    F0[i,j] <- out$par["f0"]
    F1[i,j] <- out$par["f1"]
  }

# For the selected topic (i), compare f1 estimates against the Poisson
# rates used to simulate the data.
i <- 3
e <- 1e-3
plot(dat$F[,i] + e,F1[,i] + e,pch = 20,log = "xy",
     xlab = "true f1",ylab = "estimated f1")
abline(a = 0,b = 1,col = "dodgerblue",lty = "dotted")

# For the selected topic (i), compare f0 estimates against the Poisson
# rates used to simulate the data. Since there is no "f0" used to
# simulate the data, here we approximate f0 by taking the average of
# the Poisson rates across all topics other than topic i.
plot(rowMeans(dat$F[,-i]) + e,F0[,i] + e,pch = 20,log = "xy",
     xlab = "true f0 (approx)",ylab = "estimated f0")
abline(a = 0,b = 1,col = "dodgerblue",lty = "dotted")


