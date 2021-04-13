# A short script to check that fit_poisson_glm works.

# Simulate a Poisson data set.
set.seed(1)
n  <- 200
f1 <- 0.1
f2 <- 1
s  <- sample(10,n,replace = TRUE)
q  <- runif(n)
u  <- (1-q)*f1 + q*f2
x  <- rpois(n,s*u)

# Fit the generalized linear model.
control <- glm.control(epsilon = 1e-10, maxit = 100)
L   <- cbind(s*(1-q),s*q)
dat <- data.frame(x = x,f1 = L[,1],f2 = L[,2])
fit <- glm(x ~ f1 + f2 - 1,family = poisson(link = "identity"),
           data = dat,start = c(0.5,0.5),control = control)
print(coef(fit))

# Fit the model parameters using glm with family = poisson(link =
# "identity").
out <- fit_poisson_glm(x,L)
print(out$coef)

# Compute the covariance of log(f).
print(compute_poisson_covariance(x,s*L,out$coef))
