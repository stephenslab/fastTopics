# Simulate a 80 x 100 data set.
set.seed(1)
out <- generate_test_data(80,100,3)
X   <- out$X
F   <- out$F
L   <- out$L

# Run 20 EM updates to find a good initialization.
fit0 <- fit_poisson_nmf(X,k = 3,numiter = 20,verbose = FALSE)

# COMPARE DIFFERENT METHODS
# -------------------------
fit.em     <- fit_poisson_nmf(X,fit0 = fit0,numiter = 750,method = "em",verbose = FALSE)
fit.mu     <- fit_poisson_nmf(X,fit0 = fit0,numiter = 999,method = "mu",verbose = FALSE)
fit.ccd    <- fit_poisson_nmf(X,fit0 = fit0,numiter = 300,method = "ccd",verbose = FALSE)
fit.scd    <- fit_poisson_nmf(X,fit0 = fit0,numiter = 200,method = "scd",verbose = FALSE)
fit.altsqp <- fit_poisson_nmf(X,fit0 = fit0,numiter = 150,method = "altsqp",verbose = FALSE)

# Plot improvement in solution over time.
clrs <- c("skyblue","royalblue","tomato","orange","magenta")
plot_progress_poisson_nmf(list(em = fit.em,mu = fit.mu,ccd = fit.ccd,scd = fit.scd,
                               altsqp = fit.altsqp),
                          color = clrs)

# COMPARE SPARSE VS. DENSE
# ------------------------
Y <- as(X,"dgCMatrix")
fit.em2     <- fit_poisson_nmf(Y,fit0 = fit0,numiter = 750,method = "em",verbose = FALSE)
fit.ccd2    <- fit_poisson_nmf(Y,fit0 = fit0,numiter = 300,method = "ccd",verbose = FALSE)
fit.scd2    <- fit_poisson_nmf(Y,fit0 = fit0,numiter = 200,method = "scd",verbose = FALSE)
fit.altsqp2 <- fit_poisson_nmf(Y,fit0 = fit0,numiter = 150,method = "altsqp",verbose = FALSE)

plot_progress_poisson_nmf(list(em = fit.em,mu = fit.mu,ccd = fit.ccd,scd = fit.scd,
                               altsqp = fit.altsqp,em.sp = fit.em2,ccd.sp = fit.ccd2,
                               scd.sp = fit.scd2,altsqp.sp = fit.altsqp2),
                          add.point.every = Inf,
                          color = rep(clrs,times = 2),
                          linetype = rep(c("solid","twodash"),each = 5),
                          linesize = rep(0.5,10))

# COMPARE WITH VS. WITHOUT EXTRAPOLATION
# --------------------------------------
# TO DO.

