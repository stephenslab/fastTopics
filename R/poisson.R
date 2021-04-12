# Fit the Poisson glm in which the Poisson rates are l1*f1 + ... +
# lk*fk. Inputs x and L are the data; x is a vector of length n, and L
# is an n x k matrix. Input argument f is a vector of length k giving
# initial estimates of the coefficients f1,...,fk.
#
#' @importFrom stats formula
#' @importFrom stats glm
#' @importFrom stats glm.control
#' @importFrom stats summary.glm
#' @importFrom stats poisson
#'
fit_poisson_glm <-
  function (x, L, f = rep(0.5,ncol(L)),
            control = glm.control(epsilon = 1e-10, maxit = 100)) {
  k <- length(f)
  dat <- as.data.frame(cbind(x,L))
  names(dat) <- c("x",paste0("f",1:k))
  model <- formula(paste("x ~",paste(names(dat)[-1],collapse = " + "),"- 1"))
  fit <- suppressWarnings(glm(model,family = poisson(link = "identity"),
                              data = dat,start = f,control = control))
  return(list(fit = fit,coef = summary.glm(fit)$coefficients[,"Estimate"]))
}
