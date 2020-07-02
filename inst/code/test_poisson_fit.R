# A short script to check that fit_poisson_optim and fit_poisson_em
# give the correct result.

# Simulate a Poisson data set.
set.seed(1)
n  <- 2000 # 200
f0 <- 0.1
f1 <- 1
s  <- rep(1,n) # sample(10,n,replace = TRUE)
q  <- runif(n)
u  <- (1-q)*f0 + q*f1
x  <- rpois(n,s*u)

# Fit the model parameters, f0 and f1, using optim.
out1 <- fit_poisson_optim(x,s,q)

# Fit the model parameters, f0 and f1, using the EM algorithm.
out2 <- fit_poisson_em(x,s,q)

# Fit the model parameters, f0 and f1, using the C++ implementation of
# the EM algorithm, and compare against the fit_poisson_em output.
out3 <- fit_poisson_em_rcpp(x,s,q,1,1,1e-15,40)
print(max(abs(out2$loglik - out3$loglik)))

# Fit the model parameters, f0 and f1, using the C++ implementation of
# the EM algorithm that is better suited for sparse counts.
i    <- which(x > 0)
out4 <- fit_poisson_em_sparse_rcpp(x[i],s[i],q[i],sum(s*(1-q)),sum(s*q),
                                   1,1,1e-15,40)
print(max(abs(out2$loglik - out4$loglik)))

# Compare the estimates against the values used to simulate the data.
print(data.frame(true        = c(f0,f1),
                 optim       = out1$par,
                 em          = out2$f,
                 rcpp        = with(out3,c(f0,f1)),
                 sparse_rcpp = with(out4,c(f0,f1))))

# Compare the log-likelihood at each of the solutions.
cat(sprintf("optim: %0.2f\n",-out1$value))
cat(sprintf("EM:    %0.2f\n",max(out2$loglik)))

# Compute z-score for b in ...
dat <- data.frame(x = x,q = q)
fit <- glm(x ~ q,family = poisson(link = "identity"),data = dat,
           control = list(epsilon = 1e-15,maxit = 1000))

# Compute log-fold change statistic and z-score.
# TO DO.
