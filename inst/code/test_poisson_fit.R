# A short script to check that fit_poisson_optim and fit_poisson_em
# give the correct result.

# Simulate a Poisson data set.
set.seed(1)
n  <- 200
f0 <- 0.1
f1 <- 1
s  <- sample(10,n,replace = TRUE)
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

# Finally, fit the model parameters using glm with family =
# poisson(link = "identity").
f0  <- out2$f["f0"]
f1  <- out2$f["f1"]
dat <- data.frame(x = x,b0 = s,b = s*q)
fit <- glm(x ~ b0 + b - 1,family = poisson(link = "identity"),data = dat,
           start = c(f0,f1 - f0),control = list(epsilon = 1e-15,maxit = 100))
b0  <- coef(fit)["b0"]
b   <- coef(fit)["b"]

# Compare the estimates against the values used to simulate the data.
print(data.frame(true        = c(f0,f1),
                 glm         = c(b0,b),
                 optim       = out1$par,
                 em          = out2$f,
                 rcpp        = with(out3,c(f0,f1)),
                 sparse_rcpp = with(out4,c(f0,f1)),
                 row.names = c("f0","f1")))

# Compare the log-likelihood at each of the solutions.
cat(sprintf("optim: %0.2f\n",-out1$value))
cat(sprintf("EM:    %0.2f\n",max(out2$loglik)))

# Calculate the z-scores for the glm (with identity link)
# parameterization, and compare against the internal glm calculations.
u  <- b0 + q*b
se <- sqrt(diag(solve(rbind(c(sum(x/u^2),sum(x*q/u^2)),
                            c(sum(x*q/u^2),sum(x*(q/u)^2))))))
z  <- coef(fit)/se
cat("standard errors:\n")
print(data.frame(glm = summary(fit)$coefficients[,"Std. Error"],se = se))
cat("z-scores:\n")
print(data.frame(glm = summary(fit)$coefficients[,"z value"],z = z))
            
# Calculate the z-scores for the "log-fold change" parameterization.
# TO DO.

# Plot the likelihood surface for the glm parameterization
dat <- expand.grid(list(b0 = seq(0.05,0.15,0.005),
                        b  = seq(0.8,1.1,0.005)))
dat$loglik <- 0
ns <- nrow(dat)
for (i in 1:ns) {
  b0 <- dat[i,"b0"]
  b  <- dat[i,"b"]
  dat[i,"loglik"] <- loglik_poisson(x,s*(b0 + q*b))
}
p1 <- ggplot(dat,aes(x = b0,y = b,z = loglik)) +
  geom_contour(color = "black",bins = 20) +
  theme_cowplot()

f0 <- 0.11
f1 <- 1.2

# First-order derivatives w.r.t. f0, f1.
f <- function (par) {
  u <- get_poisson_rates(s,q,par[1],par[2])
  return(loglik_poisson(x,u,1e-15))
}
u <- get_poisson_rates(s,q,f0,f1)/s
y <- (x/u - s)
c(sum(y*(1-q)),
  sum(y*q))
pracma::grad(f,c(f0,f1))

# First-order derivatives w.r.t. f0, beta.
f <- function (par) {
  u <- par[1]*(1 - q*(1 - exp(par[2])))
  return(loglik_poisson(x,s*u,1e-15))
}
beta <- log(f1/f0)
u <- f0*(1 - q*(1 - exp(beta)))
c(sum(x - s*u)/f0,
  f1*sum(q*(x/u - s)))
pracma::grad(f,c(f0,beta))

# Second-order derivatives w.r.t. f0, beta.
g <- function (par) {
  u <- par[1]*(1 - q*(1 - exp(par[2])))
  return(sum(x - s*u)/par[1])
}
c(-sum(x)/f0^2,
  -f1/f0*sum(s*q))
pracma::grad(g,c(f0,beta))

g <- function (par) {
  u  <- par[1]*(1 - q*(1 - exp(par[2])))
  f1 <- par[1]*exp(par[2])
  return(f1*sum(q*(x/u - s)))
}
c(-f1/f0*sum(s*q),
  -f1*sum(q*(s - f0*x*(1-q)/u^2)))
pracma::grad(g,c(f0,beta))

# At MLE:
-f1^2*sum(x*(q/u)^2)
