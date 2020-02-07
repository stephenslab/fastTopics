library(ggplot2)

# Simulate a 80 x 100 data set.
set.seed(1)
out <- generate_test_data(80,100,3)
X   <- out$X
F   <- out$F
L   <- out$L

# Run 20 EM updates to find a good initialization.
fit0 <- fit_poisson_nmf(X,k = 3,numiter = 50,verbose = FALSE)

# COMPARE DIFFERENT METHODS
# -------------------------
fit.mu     <- fit_poisson_nmf(X,fit0 = fit0,numiter = 500,method = "mu",verbose = FALSE)
fit.em     <- fit_poisson_nmf(X,fit0 = fit0,numiter = 400,method = "em",verbose = FALSE)
fit.ccd    <- fit_poisson_nmf(X,fit0 = fit0,numiter = 300,method = "ccd",verbose = FALSE)
fit.scd    <- fit_poisson_nmf(X,fit0 = fit0,numiter = 200,method = "scd",verbose = FALSE)
fit.altsqp <- fit_poisson_nmf(X,fit0 = fit0,numiter = 100,method = "altsqp",verbose = FALSE)

# Plot improvement in solution over time.
clrs <- c("royalblue","skyblue","tomato","orange","magenta")
p1 <- plot_progress_poisson_nmf(list(mu      = fit.mu,
                                     em     = fit.em,
                                     ccd    = fit.ccd,
                                     scd    = fit.scd,
                                     altsqp = fit.altsqp),
                                color = clrs)

# COMPARE SPARSE VS. DENSE
# ------------------------
Y <- as(X,"dgCMatrix")
fit.em2     <- fit_poisson_nmf(Y,fit0 = fit0,numiter = 500,method = "em",verbose = FALSE)
fit.ccd2    <- fit_poisson_nmf(Y,fit0 = fit0,numiter = 300,method = "ccd",verbose = FALSE)
fit.scd2    <- fit_poisson_nmf(Y,fit0 = fit0,numiter = 200,method = "scd",verbose = FALSE)
fit.altsqp2 <- fit_poisson_nmf(Y,fit0 = fit0,numiter = 150,method = "altsqp",verbose = FALSE)

p2 <- plot_progress_poisson_nmf(list(em        = fit.em,
                                     ccd       = fit.ccd,
                                     scd       = fit.scd,
                                     altsqp    = fit.altsqp,
                                     em.sp     = fit.em2,
                                     ccd.sp    = fit.ccd2,
                                     scd.sp    = fit.scd2,
                                     altsqp.sp = fit.altsqp2),
                                add.point.every = Inf,
                                color    = rep(clrs[-1],times = 2),
                                linetype = rep(c("solid","twodash"),each = 4))

# COMPARE WITH & WITHOUT EXTRAPOLATION
# ------------------------------------
fit.ccd3    <- fit_poisson_nmf(X,fit0 = fit0,numiter = 300,method = "ccd",
                               control = list(extrapolate = TRUE),verbose = FALSE)
fit.scd3    <- fit_poisson_nmf(X,fit0 = fit0,numiter = 200,method = "scd",
                               control = list(extrapolate = TRUE),verbose = FALSE)
fit.altsqp3 <- fit_poisson_nmf(X,fit0 = fit0,numiter = 100,method = "altsqp",
                               control = list(extrapolate = TRUE),verbose = FALSE)

clrs <- c("tomato","gold","magenta","firebrick","darkorange","darkmagenta")
p3 <- plot_progress_poisson_nmf(list(ccd       = fit.ccd,
                                     scd       = fit.scd,
                                     altsqp    = fit.altsqp,
                                     ccd.ex    = fit.ccd3,
                                     scd.ex    = fit.scd3,
                                     altsqp.ex = fit.altsqp3),
                                color    = clrs,
                                linesize = rep(c(0.5,0.75),each = 3))

# fit <- fit_poisson_nmf(X,fit0 = fit0,numiter = 50,method = "altsqp",
#                        control = list(extrapolate = TRUE))
