# TO DO: Explain here what this script does, and how to use it.

# SIMULATE DATA
# -------------
set.seed(1)
n <- 10
w <- rpois(n,2)
a <- abs(rnorm(n))
b <- abs(rnorm(n))

stop()

# Solve the optimization problem
#
#   minimize f(x) = sum(x*b - w*log(a + b*b))
# 
# using a simple sequential quadratic programming method.
numiter <- 100
f       <- rep(0,numiter)
x       <- 1
for (i in 1:numiter) {

  # Compute the value of the objective at x.
  y <- a + b*x;
  f <- sum(b*x - w*log(a + b*x))
    
  # Compute the gradient and Hessian at x.
  u <- l/y
  h <- sum(w*u^2)
  g <-  - sum(w*u) - u + h*x;
}
