# Verify gradient and Hessian calculations for the "single gene"
# Poisson model.
library(pracma)

# Simulate data x ~ Pois(u), with u = l0*f0 + l1*f1.
set.seed(1)
n  <- 40
f0 <- 0.1
f1 <- 1
s  <- sample(10,n,replace = TRUE)
u  <- runif(n)
l0 <- s*(1-u)
l1 <- s*u
x  <- rpois(n,l0*f0 + l1*f1)

# Compute the log-likelihood under the model x ~ Pois(u), with
# Poisson rates u = l0*f0 + l1*f1.
loglik <- function (x, l0, l1, f0, f1)
  sum(dpois(x,l0*f0 + l1*f1,log = TRUE))

# Compute the gradient of the log-likelihood with respect to log(f0)
# and log(f1).
loglik_grad <- function (x, l0, l1, f0, f1) {
  u <- l0*f0 + l1*f1
  y <- x/u - 1
  return(c(f0*sum(l0*y),
           f1*sum(l1*y)))
}

# Compute the 2 x 2 matrix of second-order partial derivatives of the
# log-likelihood with respect to log(f0) and log(f2).
loglik_hessian <- function (x, l0, l1, f0, f1) {
  u <- l0*f0 + l1*f1
  y <- x/u - 1
  z <- x/u^2
  return(rbind(c(f0*sum(l0*y) - f0^2*sum(l0^2*z),-f0*f1*sum(l0*l1*z)),
               c(-f0*f1*sum(l0*l1*z),f1*sum(l1*y) - f1^2*sum(l1^2*z))))
}

# Compare loglik_grad against numerical gradients calculated using
# finite differences.
f0 <- runif(1)
f1 <- runif(1)

cat("gradient:\n")
print(grad(function (v) loglik(x,l0,l1,exp(v[1]),exp(v[2])),log(c(f0,f1))),
      digits = 12)
print(loglik_grad(x,l0,l1,f0,f1),digits = 12)

cat("Hessian:\n")
print(rbind(grad(function (v) loglik_grad(x,l0,l1,exp(v[1]),exp(v[2]))[1],
                 log(c(f0,f1))),
            grad(function (v) loglik_grad(x,l0,l1,exp(v[1]),exp(v[2]))[2],
                 log(c(f0,f1)))),
      digits = 12)
print(loglik_hessian(x,l0,l1,f0,f1),digits = 12)
