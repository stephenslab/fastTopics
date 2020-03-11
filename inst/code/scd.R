# A small script to illustrate the co-ordinate ascent updates in the
# sequential co-ordinate descent (SCD) algorithm described by Lin &
# Boutros (2018).

# SIMULATE DATA
# -------------
set.seed(49)
n <- 20
w <- rpois(n,2)
a <- abs(rnorm(n))
b <- abs(rnorm(n))

# Solve the following 1-d optimization problem:
#
#   minimize f(x) = sum(b*x - w*log(y))
#   subject to y = a + b*x,
#              x >= 0.
# 
# using a simple sequential quadratic programming (SQP) method.
x       <- 1
numiter <- 20
f       <- rep(0,numiter)
for (i in 1:numiter) {

  # Compute the value of the objective at x.
  y    <- a + b*x;
  f[i] <- sum(b*x - w*log(y))
    
  # Compute the gradient and Hessian at x.
  u <- b/y
  h <- sum(w*u^2)
  g <- sum(b - w*u)

  # Update x.
  x <- max(0,x - g/h)
}
cat(sprintf("solution: %0.6f\n",x))

# Plot the improvement in the solution over time.
plot(1:numiter,f - min(f) + 1e-15,type = "b",pch = 20,col = "darkblue",
     log = "y",xlab = "iteration",ylab = "distance to minimum")

