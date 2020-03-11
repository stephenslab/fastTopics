# A small script to illustrate the co-ordinate ascent updates in the
# sequential co-ordinate descent (SCD) algorithm described by Lin &
# Boutros (2018).
#
# Here I enhance the SCD algorithm with a simple backtracking line
# search to guarantee that the objective decreases at each iteration.
#

# SCRIPT PARAMETERS
# -----------------
n           <- 20
numiter     <- 20
line.search <- TRUE

# SIMULATE DATA
# -------------
set.seed(49)
w <- rpois(n,2)
a <- abs(rnorm(n))
b <- abs(rnorm(n))

# Solve the following 1-d optimization problem:
#
#   minimize    f(x) = sum(b*x - w*log(y))
#   subject to  y = a + b*x,
#               x >= 0.
# 
# using a simple sequential quadratic programming (SQP) method.
x <- 1
e <- 1e-15
f <- rep(0,numiter)
for (i in 1:numiter) {

  # Compute the value of the objective at x.
  y    <- a + b*x;
  f[i] <- sum(b*x - w*log(y))
    
  # Compute the gradient and Hessian at x.
  u <- b/y
  h <- sum(w*u^2)
  g <- sum(b - w*u)

  # Optionally, perform backtracking line search to determine a
  # suitable step size.
    p <- -g/h
  if (line.search) {
    if (p >= -e)
      s <- 1
    else
      s <- min(1,-x/p)
    smin <- e
    while (TRUE) {
      xnew <- x + s*p
      ynew <- a + b*xnew
      fnew <- sum(b*xnew - w*log(ynew))
      if (s < smin) {
        xnew <- x
        s    <- 0
        break
      } else if (fnew < f[i])
        break
      else
        s <- s/2
    }
  } else
    xnew <- max(0,x + p)
      
  # Update x.
  x <- xnew
}
cat(sprintf("solution: %0.6f\n",x))

# Plot the improvement in the solution over time.
y <- f - min(f) + 1e-15
plot(1:numiter,y,type = "l",col = "dodgerblue",lwd = 1,log = "y",
     xlab = "iteration",ylab = "distance to minimum")
points(1:numiter,y,pch = 20,col = "dodgerblue")

# This is a plot I created to compare the progression with and without
# the backtracking line search.
# 
#   f <- min(c(f.nols,f.ls))
#   y <- f.nols - f + 1e-8
#   plot(1:numiter,y,type = "l",col = "dodgerblue",lwd = 1,log = "y",
#        xlab = "iteration",ylab = "distance to minimum")
#   points(1:numiter,y,pch = 20,col = "dodgerblue")
#   y <- f.ls - f + 1e-8
#   lines(1:numiter,y,col = "darkorange",lwd = 1)
#   points(1:numiter,y,pch = 20,col = "darkorange")
#
