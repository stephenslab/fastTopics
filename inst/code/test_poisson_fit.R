# Verify the Poisson model computations on a small example with k = 2.
library(ggplot2)
library(cowplot)

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
print(log(coef(fit)))

# Fit the model parameters using glm with family = poisson(link =
# "identity").
out <- fit_poisson_glm(x,L)
print(log(out$coef))

# Compute the covariance of log(f).
cat("Cov(log(f)) estimated via Laplace approximation:\n")
print(compute_poisson_covariance(x,L,out$coef))

# Draw samples from the posterior using random-walk Metropolis.
ns <- 1e5
set.seed(1)
print(system.time(
  sim <- simulate_posterior_poisson(x,L,out$coef,ns = ns,s = 0.3)))
set.seed(1)
D <- matrix(rnorm(2*ns),ns,2)
U <- matrix(runif(2*ns),ns,2)
print(system.time(
  sim2 <- simulate_posterior_poisson_rcpp(x,L,out$coef,D,U,0.3,1e-15)))

cat("MCMC estimate of Cov(log(f)):\n")
print(cov(log(sim$samples)))
cat(sprintf("Acceptance rate: %0.2f%%\n",100*sim$ar))

# Get 90% HPD intervals.
cat("MCMC estimates of 90% HPD intervals:\n")
print(hpd(log(sim$samples[,1]),0.9))
print(hpd(log(sim$samples[,2]),0.9))

# Plot the likelihood surface.
dat     <- expand.grid(t1 = seq(-4,1,0.05),t2 = seq(-4,1,0.02))
n       <- nrow(dat)
dat$lik <- 0
for (i in 1:n) {
  f <- exp(c(dat[i,1],dat[i,2]))
  u <- drop(L %*% f)
  dat[i,"lik"] <- sum(dpois(x,u,log = TRUE))
}
dat$lik <- exp(dat$lik - max(dat$lik))
p1 <- ggplot(dat,aes(x = t1,y = t2,z = lik)) +
  geom_contour(color = "black",bins = 16) +
  geom_point(data = as.data.frame(t(log(out$coef))),
             mapping = aes(x = f1,y = f2),
             color = "red",shape = 4,
             inherit.aes = FALSE) +
  labs(x = "log(f1)",y = "log(f2)") +
  theme_cowplot(font_size = 10)

# Plot the MCMC density estimate.
sim$samples <- as.data.frame(log(sim$samples))
names(sim$samples) <- c("k1","k2")
p2 <- ggplot(sim$samples,aes(x = k1,y = k2)) +
  geom_density_2d(color = "black") +
  geom_point(data = as.data.frame(t(log(out$coef))),
             mapping = aes(x = f1,y = f2),
             color = "red",shape = 4,
             inherit.aes = FALSE) +
  labs(x = "log(f1)",y = "log(f2)") +
  theme_cowplot(font_size = 10)
print(plot_grid(p1,p2))
