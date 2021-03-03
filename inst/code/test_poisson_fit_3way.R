# A short script to check that my manual calculations give the correct
# result for the Poisson mixture model with three mixture components.

# Simulate a Poisson data set.
set.seed(1)
n  <- 40
f1 <- 0.1
f2 <- 0.2
f3 <- 1
s  <- sample(10,n,replace = TRUE)
q1 <- runif(n)
q2 <- runif(n)
q3 <- runif(n)
z  <- q1 + q2 + q3
q1 <- q1/z
q2 <- q2/z
q3 <- q3/z
u  <- q1*f1 + q2*f2 + q3*f3
x  <- rpois(n,s*u)

# Fit the model using glm with family = poisson(link = "identity").
dat <- data.frame(x = x,b1 = s*q1,f2 = s*(q1 + q2),f3 = s*(q1 + q3))
fit <- suppressWarnings(glm(x ~ b1 + f2 + f3 - 1,
                            family = poisson(link = "identity"),
                            data = dat,start = c(f1 - f2 - f3,f2,f3),
                            control = glm.control(epsilon=1e-10,maxit=100)))
ans <- summary.glm(fit)$coefficients
b1  <- ans["b1","Estimate"]
f2  <- ans["f2","Estimate"]
f3  <- ans["f3","Estimate"]

stop()

  # Output the parameter estimates, log-fold changes and test statistics.
  ans <- summary.glm(fit)$coefficients
  b0  <- ans["b0","Estimate"]
  b   <- ans["b","Estimate"]
  f0  <- b0
  f1  <- b + b0
  return(c(f0   = f0,
           f1   = f1,
           beta = log2(f1/f0),
           se   = ans["b","Std. Error"],
           z    = ans["b","z value"],
           pval = -log10(ans["b","Pr(>|z|)"])))
fit <- fit_poisson_glm(x,s,q)

# Fit the model parameters, f1, f2 and f3, using optim.
out <- fit_poisson_optim(x,s,q)
