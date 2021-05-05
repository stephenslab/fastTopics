# Illustrate use of random-walk Metropolis to simulate a Gamma random
# variable.
n <- 1e5
a <- 2
s <- 4
set.seed(1)
X <- rep(0,n)
Y <- rgamma(n,shape = a,scale = s)
x <- 0
r <- 0
for (i in 1:n) {
  xnew <- x + rnorm(1,sd = 1)
  p <- min(1,exp(xnew - x +
                 dgamma(exp(xnew),shape = a,scale = s,log = TRUE) -
                 dgamma(exp(x),shape = a,scale = s,log = TRUE)))
  if (runif(1) < p) {
    x <- xnew
    r <- r + 1
  }
  X[i] <- exp(x)
}
print(mean(X))
print(mean(Y))
print(sd(X))
print(sd(Y))
print(r/n)
