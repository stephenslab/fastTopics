# TO DO: Explain here what this script is for, and how to use it.

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
cat(coef(fit),f0,"\n")

# Compute the s.e. of log(f0) in the basic Poisson model.
se <- 1/sqrt(f0*sum(s))

# Compute the s.e. of log(f0) using numerical integration.
ns <- 1000
t  <- seq(-5,-1,length.out = ns)
w  <- rep(0,ns)
for (i in 1:ns) {
  w[i] <- sum(dpois(x,s*exp(t[i]),log = TRUE))
}
w <- exp(w - max(w))
w <- w/sum(w)
mu <- sum(w*t)
se_mc <- sqrt(sum(w*(t - mu)^2))

# The two s.e. calculations should be nearly the same.
cat(se,se_mc,"\n")
