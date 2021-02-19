# A short script to check that fit_poisson_optim and fit_poisson_em
# give the correct result, and to verify my standard error (s.e.) and
# z-score calculations for the Poisson model.

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
out3 <- fit_poisson_em_rcpp(x,s,q,1,1,1e-15,100,1e-8)

# Fit the model parameters, f0 and f1, using the C++ implementation of
# the EM algorithm that is better suited for sparse counts.
i    <- which(x > 0)
out4 <- fit_poisson_em_sparse_rcpp(x[i],s[i],q[i],sum(s*(1-q)),sum(s*q),
                                   1,1,1e-15,100,1e-8)

# Compare the log-likelihood at each of the solutions.
cat(sprintf("optim:       %0.12f\n",-out1$value))
cat(sprintf("EM:          %0.12f\n",max(out2$loglik)))
cat(sprintf("EM (rcpp):   %0.12f\n",max(out3$loglik)))
cat(sprintf("EM (sparse): %0.12f\n",max(out4$loglik)))

# Finally, fit the model parameters using glm with family =
# poisson(link = "identity"). Note that the parameterization is
# slightly different: b0 = f0 and b = f1 - f0.
fit <- fit_poisson_glm(x,s,q)

# Compare the estimates against the values used to simulate the data.
print(data.frame(glm         = fit$f,
                 optim       = out1$par,
                 em          = out2$f,
                 rcpp        = with(out3,c(f0,f1)),
                 sparse_rcpp = with(out4,c(f0,f1)),
                 row.names = c("f0","f1")))

# Z-SCORE CALCULATIONS
# --------------------
# Manually calculate the z-scores for the glm (with identity link)
# parameterization, and compare against the internal glm calculations.
# They should be close.
b0   <- fit$f[1]
b    <- fit$f[2] - b0
u    <- b0 + q*b
se   <- sqrt(diag(solve(rbind(c(sum(x/u^2),sum(x*q/u^2)),
                              c(sum(x*q/u^2),sum(x*(q/u)^2))))))
se   <- se[2]
z    <- b/se
pval <- 2*pnorm(-abs(z))
cat("standard errors:\n")
print(data.frame(glm = fit$se,se = se))
cat("z-scores:\n")
print(data.frame(glm = fit$z,z = z))
cat("p-values:\n")
print(data.frame(glm = fit$pval,pval = pval))
