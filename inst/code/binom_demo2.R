# TO DO: Explain here what this script is for, and how to use it.
set.seed(1)

# Simulate data.
n   <- 80
m   <- 100
k   <- 5
dat <- simulate_count_data(n,m,k)
X   <- dat$X

# Fit multinomial topic model.
F   <- rand(m,k)
L   <- dat$L
fit <- fit_poisson_nmf(X,fit0 = init_poisson_nmf(X,F,L),numiter = 20,
                       update.loadings = NULL,verbose = FALSE)
fit <- poisson2multinom(fit)

# Fit a binomial topic model for each j and topic k.
P0opt <- matrix(0,m,k)
P1opt <- matrix(0,m,k)
P0em  <- matrix(0,m,k)
P1em  <- matrix(0,m,k)
y     <- rowSums(X)
for (i in 1:m)
  for (j in 1:k) {
    out1 <- fit_binom_optim(X[,i],y - X[,i],fit$L[,j])
    out2 <- fit_binom_optim(X[,i],y - X[,i],fit$L[,j])
    P0opt[i,j] <- out1$par["p0"]
    P1opt[i,j] <- out1$par["p1"]
    P0em[i,j]  <- out2$p["p0"]
    P1em[i,j]  <- out2$p["p1"]
  }

# Compare the optim and EM estimates of the binomial topic model
# parameters.
plot(P0opt,P0em,pch = 20)
plot(P1opt,P1em,pch = 20)

stop()

k <- 4
j <- 15
x <- X[,j]
y <- rowSums(X[,-j])
q <- fit$L[,k]
out1 <- fit_binom_optim(x,y,q)
out2 <- fit_binom_em(x,y,q,numiter = 40)
print(out1$par)

out <- compute_lfoldchange_helper(X,fit$F,fit$L,k)
n0 <- out$n0
n1 <- out$n1
p0 <- (sum(x) - n1[j])/sum(n0[-j])
p1 <- n1[j]/(sum(y) - sum(n0[-j]))
p0 <- p0/(1 + p0)
p1 <- p1/(1 + p1)


