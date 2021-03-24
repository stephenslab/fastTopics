# TO DO: Explain here what this script is for, and how to use it.
library(pracma)

# Simulate binomial data, x ~ binom(s*p0).
set.seed(1)
n  <- 50
s  <- ceiling(10*runif(n))
p0 <- 0.05
q  <- runif(n)
x  <- rbinom(n,s,p0)

# Fit the basic Poisson model x ~ Pois(s*f0) using glm.
fit <- glm(x ~ f0 - 1,family = poisson(link = "identity"),
           data = data.frame(x = x,f0 = s),start = 0.5,
           control = glm.control(epsilon = 1e-10, maxit = 100))

# Compute the MLE of p0 in the basic Poisson model.
f0 <- sum(x)/sum(s)

# The glm estimate should be the same as f0.
print(coef(fit) - f0)

# Compute the s.e. of log(f0) in the basic Poisson model.
# TO DO.

# Compute the s.e. of log(f0) using numerical integration.
# TO DO.
