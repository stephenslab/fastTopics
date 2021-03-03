# A short script to check that my manual calculations give the correct
# result for the Poisson mixture model with three mixture components.

# Simulate a Poisson data set.
set.seed(1)
n  <- 40
f1 <- 0.1
f2 <- 0.2
f3 <- 1
s  <- sample(10,n,replace = TRUE)
q1 <- runif(n,min = 0,max = 1)
q2 <- runif(n,min = 0,max = 0.5)
q3 <- runif(n,min = 0,max = 0.1)
z  <- q1 + q2 + q3
q1 <- q1/z
q2 <- q2/z
q3 <- q3/z
u  <- q1*f1 + q2*f2 + q3*f3
x  <- rpois(n,s*u)

# Fit the reparameterized model using glm with family = poisson(link =
# "identity").
dat <- data.frame(x = x,b1 = s*q1,f2 = s*(q1 + q2),f3 = s*(q1 + q3))
fit <- suppressWarnings(glm(x ~ b1 + f2 + f3 - 1,
                            family = poisson(link = "identity"),
                            data = dat,start = c(f1 - f2 - f3,f2,f3),
                            control = glm.control(epsilon=1e-10,maxit=100)))
ans <- summary.glm(fit)$coefficients
b1  <- ans["b1","Estimate"]
f2  <- ans["f2","Estimate"]
f3  <- ans["f3","Estimate"]

# Fit the model using optim.
e <- 1e-15
loglik_poisson <- function (x, y)
  return(sum(x*log(y + e) - y))
get_poisson_rates <- function (q1, q2, q3, f1, f2, f3)
  q1*f1 + q2*f2 + q3*f3
f <- function (par) {
  u <- get_poisson_rates(q1,q2,q3,par[1],par[2],par[3])
  return(-loglik_poisson(x,s*u))
}
g <- function (par) {
  u <- get_poisson_rates(q1,q2,q3,par[1],par[2],par[3])
  y <- (s*u - x)/(u + e)
  return(c(sum(y*q1),sum(y*q2),sum(y*q3)))
}
out <- optim(c(f1,f2,f3),f,g,method = "L-BFGS-B",lower = c(e,e,e),
             control = list(factr = 1e5, maxit = 100))
out  <- within(out,{
                 f1 <- par[1]
                 f2 <- par[2]
                 f3 <- par[3]
                 b1 <- f1 - f2 - f3
               })

# Compare the glm and optim parameter estimates.
print(matrix(c(b1,f2,f3,
               out$b1,out$f2,out$f3),2,3,byrow = TRUE,
             dimnames = list(c("glm","optim"),c("b1","f2","f3"))),
      digits = 8)

# Manually calculate the standard errors, z-scores and p-values.
f1  <- out$f1
f2  <- out$f2
f3  <- out$f3
u   <- get_poisson_rates(q1,q2,q3,f1,f2,f3)
h11 <- sum(x*(q1/u)^2)
h22 <- sum(x*(q2/u)^2)
h33 <- sum(x*(q3/u)^2)
h12 <- sum(x*q1*q2/u^2)
h13 <- sum(x*q1*q3/u^2)
h23 <- sum(x*q2*q3/u^2)
B <- rbind(c(1,-1,-1),
           c(0, 1, 0),
           c(0, 0, 1))
H <- rbind(c(h11,h12,h13),
           c(h12,h22,h23),
           c(h13,h23,h33))
se   <- sqrt(diag(B %*% solve(H) %*% t(B)))
z    <- c(f1 - f2 - f3,f2,f3)/se
pval <- 2*pnorm(-abs(z))

# Compare the glm statistics against my manual calculations.
cat("standard errors:\n")
print(cbind(glm = ans[,"Std. Error"],optim = se))
cat("z-scores:\n\n")
print(cbind(glm = ans[,"z value"],optim = z))
cat("p-values:\n")
print(cbind(glm = ans[,"Pr(>|z|)"],optim = pval))
