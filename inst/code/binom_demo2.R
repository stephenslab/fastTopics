# TO DO: Explain here what this script is for, and how to use it.
source("binom.R")
set.seed(1)

# Simulate data.
n   <- 200
m   <- 400
k   <- 3
dat <- simulate_count_data(n,m,k)
X   <- dat$X

# Fit multinomial topic model.
F   <- rand(m,k)
L   <- dat$L
fit <- fit_poisson_nmf(X,fit0 = init_poisson_nmf(X,F,L),numiter = 20,
                       update.loadings = NULL,verbose = FALSE)

# Fit a binomial topic model for each j and topic k.
fit.optim <- fit_binom_topic_model(X,poisson2multinom(fit),method = "optim")
# fit.em <- fit_binom_topic_model(...)

stop()

# Compute the MLEs for the binomial topic model directly from the MLEs
# of the multinomial topic model.
out <- multinom2binom(X,fit,e = 0,version = "R")

# Compare the optim and EM estimates of the binomial topic model
# parameters.
plot(P0opt,P0em,pch = 20)
plot(P1opt,P1em,pch = 20)

# Compare the optim estimates to the multinom2binom calculations.
plot(P0opt,out$P0,pch = 20)
abline(a = 0,b = 1,col = "magenta",lty = "dotted")
plot(P1opt,out$P1,pch = 20)
abline(a = 0,b = 1,col = "magenta",lty = "dotted")
