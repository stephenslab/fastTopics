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

# Fit a Poisson model to each combination of gene j and topic k.
F0optim      <- matrix(0,m,k)
F1optim      <- matrix(0,m,k)
F0em         <- matrix(0,m,k)
F1em         <- matrix(0,m,k)
loglik.optim <- matrix(0,m,k)
loglik.em    <- matrix(0,m,k)
for (i in 1:m)
  for (j in 1:k) {
    out1              <- fit_poisson_optim(X[,i],s,L[,j])
    out2              <- fit_poisson_em(X[,i],s,L[,j])
    F0optim[i,j]      <- out1$par["f0"]
    F1optim[i,j]      <- out1$par["f1"]
    F0em[i,j]         <- out2$f["f0"]
    F1em[i,j]         <- out2$f["f1"]
    loglik.optim[i,j] <- -out1$value
    loglik.em[i,j]    <- max(out2$loglik)
  }

# Check that EM and optim produce the same, or nearly the same,
# estimates of the model parameters, f0 and f1.
e <- 1e-3
plot(F1optim + e ,F1em + e,pch = 20,log = "xy")
abline(a = 0,b = 1,col = "dodgerblue",lty = "dotted")

# Compare the EM and optim log-likelihoods.
print(quantile(loglik.em - loglik.optim,seq(0,1,0.2)))

# For the selected topic (i), compare f1 estimates against the Poisson
# rates used to simulate the data.
i <- 3
plot(dat$F[,i] + e,F1em[,i] + e,pch = 20,log = "xy",
     xlab = "true f1",ylab = "estimated f1")
abline(a = 0,b = 1,col = "dodgerblue",lty = "dotted")

# For the selected topic (i), compare f0 estimates against the Poisson
# rates used to simulate the data. Since there is no "f0" used to
# simulate the data, here we approximate f0 by taking the average of
# the Poisson rates across all topics other than topic i.
plot(rowMeans(dat$F[,-i]) + e,F0em[,i] + e,pch = 20,log = "xy",
     xlab = "true f0 (approx)",ylab = "estimated f0")
abline(a = 0,b = 1,col = "dodgerblue",lty = "dotted")


