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

# Compare the estimates against the values used to simulate the data.
print(data.frame(true        = c(f0,f1),
                 glm         = c(b0,b),
                 optim       = out1$par,
                 em          = out2$f,
                 rcpp        = with(out3,c(f0,f1)),
                 sparse_rcpp = with(out4,c(f0,f1))))

# Compare the log-likelihood at each of the solutions.
cat(sprintf("optim: %0.2f\n",-out1$value))
cat(sprintf("EM:    %0.2f\n",max(out2$loglik)))

# Calculate the z-scores for the glm (with identity link)
# parameterization, and compare against the internal glm calculations.
b0 <- coef(fit)["b0"]
b  <- coef(fit)["b"]
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
