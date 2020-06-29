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
t1 <- system.time(out.em  <- fit_univar_poisson_models(X,L,s,method="em"))
t2 <- system.time(out.em2 <- fit_univar_poisson_models(X,L,s,method="em-rcpp"))
t3 <- system.time(out.optim <- fit_univar_poisson_models(X,L,s,method="optim"))
print(t1)
print(t2)
print(t3)

# Check that the R and C++ implementations of the EM algorithm produce
# the same, or nearly the same, likelihoods and estimates of the model
# parameters, f0 and f1.
print(max(abs(out.em$F0 - out.em2$F0)))
print(max(abs(out.em$F1 - out.em2$F1)))
print(max(abs(out.em$loglik - out.em2$loglik)))

# Check that EM and optim produce the same, or nearly the same,
# estimates of the model parameters, f0 and f1.
e <- 1e-3
plot(out.optim$F1 + e,out.em$F1 + e,pch = 20,log = "xy")
abline(a = 0,b = 1,col = "dodgerblue",lty = "dotted")

# Compare the EM and optim log-likelihoods.
print(range((out.em$loglik - out.optim$loglik)/out.optim$loglik))

# For the selected topic (i), compare f1 estimates against the Poisson
# rates used to simulate the data.
i <- 3
plot(dat$F[,i] + e,out.em$F1[,i] + e,pch = 20,log = "xy",
     xlab = "true f1",ylab = "estimated f1")
abline(a = 0,b = 1,col = "dodgerblue",lty = "dotted")

# For the selected topic (i), compare f0 estimates against the Poisson
# rates used to simulate the data. Since there is no "f0" used to
# simulate the data, here we approximate f0 by taking the average of
# the Poisson rates across all topics other than topic i.
plot(rowMeans(dat$F[,-i]) + e,out.em$F0[,i] + e,pch = 20,log = "xy",
     xlab = "true f0 (approx)",ylab = "estimated f0")
abline(a = 0,b = 1,col = "dodgerblue",lty = "dotted")


