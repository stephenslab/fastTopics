#' @rdname predict
#'
#' @title Title Goes Here
#'
#' @description Description goes here.
#'
#' @param object Describe input argument "object" here.
#'
#' @param newdata Describe (optional) input "newdata" here.
#'
#' @param numiter Describe input "numiter" here.
#'
#' @param \dots Describe ... here.
#'
#' @seealso \code{\link{fit_poisson_nmf}}
#' 
#' @examples
#' # Simulate a 175 x 1,200 counts matrix.
#' set.seed(1)
#' dat <- simulate_count_data(175,1200,k = 3)
#'
#' # Split the data into training and test sets.
#' train <- dat$X[1:100,]
#' test <- dat$X[101:175,]
#'
#' # Fit a Poisson non-negative matrix factorization using the
#' # training data.
#' fit <- init_poisson_nmf(train,F = dat$F,init.method = "random")
#' fit <- fit_poisson_nmf(train,fit0 = fit)
#'
#' # Compare the estimated loadings in the training data against the
#' # loadings used to simulate these data.
#' Ltrain <- predict(fit)
#' plot(dat$L[1:100,],Ltrain,pch = 20,col = "darkblue")
#' abline(a = 0,b = 1,col = "magenta",lty = "dotted",
#'        xlab = "true",ylab = "estimated")
#'
#' # Next, predict loadings in unseen (test) data points, and compare
#' # these predictions against the loadings that were used to simulate
#' # the test data.
#' Ltest <- predict(fit,test)
#' plot(dat$L[101:175,],Ltest,pch = 20,col = "darkblue",
#'      xlab = "true",ylab = "estimated")
#' abline(a = 0,b = 1,col = "magenta",lty = "dotted")
#'
#' # Simulate a 175 x 1,200 counts matrix.
#' set.seed(1)
#' dat <- simulate_multinom_gene_data(175,1200,3)
#'
#' # Split the data into training and test sets.
#' train <- dat$X[1:100,]
#' test <- dat$X[101:175,]
#'
#' # Fit a topic model using the training data.
#' fit <- init_poisson_nmf(train,F = dat$F,init.method = "random")
#' fit <- fit_poisson_nmf(train,fit0 = fit)
#' fit <- poisson2multinom(fit)
#'
#' # Compare the estimated topic proportions in the training data against
#' # the topic proportions used to simulate these data.
#' Ltrain <- predict(fit)
#' plot(dat$L[1:100,],Ltrain,pch = 20,col = "darkblue")
#' abline(a = 0,b = 1,col = "magenta",lty = "dotted",
#'        xlab = "true",ylab = "estimated")
#'
#' # Next, predict loadings in unseen (test) data points, and compare
#' # these predictions against the loadings that were used to simulate
#' # the test data.
#' Ltest <- predict(fit,test)
#' plot(dat$L[101:175,],Ltest,pch = 20,col = "darkblue",
#'      xlab = "true",ylab = "estimated")
#' abline(a = 0,b = 1,col = "magenta",lty = "dotted")
#'
#' @importFrom stats predict
#' 
#' @method predict poisson_nmf_fit
#'
#' @export
#' 
predict.poisson_nmf_fit <- function (object, newdata, numiter = 20, ...) {
  if (missing(newdata))
    return(object$L)
  else
    return(project_poisson_nmf(object,newdata,numiter,...)$L)
}

#' @rdname predict
#' 
#' @importFrom stats predict
#' 
#' @method predict multinom_topic_model_fit
#'
#' @export
#' 
predict.multinom_topic_model_fit <- function (object, newdata, numiter = 20,
                                              ...) {
  if (missing(newdata))
    return(object$L)
  else
    return(poisson2multinom(project_poisson_nmf(multinom2poisson(object),
                                                newdata,numiter,...))$L)
}

# TO DO: Explain here what this function does, and how to use it.
project_poisson_nmf <- function (fit0, X, numiter, ...) {
  n   <- nrow(X)
  k   <- ncol(fit0$F)
  L0  <- matrix(1/k,n,k)
  fit <- init_poisson_nmf(X,F = fit0$F,L = L0)
  fit <- fit_poisson_nmf(X,fit0 = fit,numiter = numiter,...)
  return(fit)
}
