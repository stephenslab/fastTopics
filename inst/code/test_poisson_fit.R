# A short script to check that fit_poisson_optim and fit_poisson_em
# give the correct result, and to verify my standard error (s.e.) and
# z-score calculations for the Poisson model.
library(pracma)

# Simulate a Poisson data set.
set.seed(1)
n  <- 40
f0 <- 0.1
f1 <- 1
s  <- sample(10,n,replace = TRUE)
q  <- runif(n)
u  <- (1-q)*f0 + q*f1
x  <- rpois(n,s*u)

# COMPUTE MLEs OF f0,f1
# ---------------------
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
# poisson(link = "identity"). Note that the parameterization is
# slightly different: b0 = f0 and b = f1 - f0.
dat <- data.frame(x = x,b0 = s,b = s*q)
fit <- glm(x ~ b0 + b - 1,family = poisson(link = "identity"),data = dat,
           start = c(f0,f1 - f0))
b0  <- coef(fit)["b0"]
b   <- coef(fit)["b"]

# Compare the estimates against the values used to simulate the data.
print(data.frame(glm         = c(b0,b + b0),
                 optim       = out1$par,
                 em          = out2$f,
                 rcpp        = with(out3,c(f0,f1)),
                 sparse_rcpp = with(out4,c(f0,f1)),
                 row.names = c("f0","f1")))

# Compare the log-likelihood at each of the solutions.
cat(sprintf("optim:       %0.6f\n",-out1$value))
cat(sprintf("EM:          %0.6f\n",max(out2$loglik)))
cat(sprintf("EM (rcpp):   %0.6f\n",max(out3$loglik)))
cat(sprintf("EM (sparse): %0.6f\n",max(out4$loglik)))

# Z-SCORE CALCULATIONS
# --------------------
# Manually calculate the z-scores for the glm (with identity link)
# parameterization, and compare against the internal glm
# calculations. They should be very similar.
u    <- b0 + q*b
se   <- sqrt(diag(solve(rbind(c(sum(x/u^2),sum(x*q/u^2)),
                              c(sum(x*q/u^2),sum(x*(q/u)^2))))))
z    <- coef(fit)/se
pval <- 2*pnorm(-abs(z))
cat("standard errors:\n")
print(data.frame(glm = summary(fit)$coefficients[,"Std. Error"],se = se))
cat("z-scores:\n")
print(data.frame(glm = summary(fit)$coefficients[,"z value"],z = z))
cat("p-values:\n")
print(data.frame(glm = summary(fit)$coefficients[,"Pr(>|z|)"],pval = pval))

# Calculate the z-scores for the "log-fold change" parameterization.
f0 <- out2$f["f0"]
f1 <- out2$f["f1"]
b  <- log(f1/f0)
u  <- get_poisson_rates(q,f0,f1)
se <- sqrt(diag(solve(rbind(c(sum(x)/f0^2,f1/f0*sum(s*q)),
                            c(f1/f0*sum(s*q),f1^2*sum(x*(q/u)^2))))))
z  <- c(f0,b)/se
pval <- 2*pnorm(-abs(z))
print(data.frame(mle  = c(f0,b/log(2)),
                 se   = se/c(1,log(2)),
                 z    = z,
                 pval = pval,
                 row.names = c("f0","beta")))
print(compute_poisson_zscore(x,q,s,f0,f1))
print(log2(f1/f0))

# GRADIENT & HESSIAN CALCULATIONS
# -------------------------------
f0 <- 0.11
f1 <- 1.2

# First-order derivatives w.r.t. f0, f1.
f <- function (par) {
  f0 <- par[1]
  f1 <- par[2]
  u  <- get_poisson_rates(q,f0,f1)
  return(loglik_poisson(x,s*u,1e-15))
}
u <- get_poisson_rates(q,f0,f1)
y <- (x/u - s)
print(cbind(c(sum(y*(1-q)),sum(y*q)),
            grad(f,c(f0,f1))))

# First-order derivatives w.r.t. f0, beta.
get_poisson_rates2 <- function (q, f0, b)
  f0*(1 - q*(1 - exp(b)))
f <- function (par) {
  f0 <- par[1]
  b  <- par[2]
  u  <- get_poisson_rates2(q,f0,b)
  return(loglik_poisson(x,s*u,1e-15))
}
b <- log(f1/f0)
u <- get_poisson_rates2(q,f0,b)
print(cbind(c(sum(x - s*u)/f0,f1*sum(q*(x/u - s))),
            grad(f,c(f0,b))))

f0 <- out2$f["f0"]
f1 <- out2$f["f1"]
b  <- log(f1/f0)

# Second-order derivatives w.r.t. f0, beta.
g <- function (par) {
  f0 <- par[1]
  b  <- par[2]
  u  <- get_poisson_rates2(q,f0,b)
  return(sum(x - s*u)/f0)
}
print(cbind(c(-sum(x)/f0^2,-f1/f0*sum(s*q)),
            grad(g,c(f0,b))))

# Second-order derivatives w.r.t. f0, beta.
g <- function (par) {
  f0 <- par[1]
  b  <- par[2]
  u  <- get_poisson_rates2(q,f0,b)
  f1 <- f0*exp(b)
  return(f1*sum(q*(x/u - s)))
}
u <- get_poisson_rates2(q,f0,b)
print(cbind(c(-f1/f0*sum(s*q),-f1*sum(q*(s - f0*x*(1-q)/u^2))),
            grad(g,c(f0,b))))

# At the MLE, the second-order derivatives simplify.
print(-f1^2*sum(x*(q/u)^2))
