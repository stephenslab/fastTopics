#' @rdname predict
#'
#' @title Predict Methods for Poisson NMF and Multinomial Topic Model
#'
#' @description Predict loadings based on previously fit Poisson NMF,
#'   or predict topic proportions based on previously fit multinomial
#'   topic model. This can be thought of as projecting data points onto
#'   a previously estimated set of factors \code{fit$F}.
#'
#' @param object An object of class \dQuote{poisson_nmf_fit} or
#'   \dQuote{multinom_topic_model_fit}.
#' 
#' @param newdata An optional counts matrix. If omitted, the loadings
#'   estimated in the original data are returned.
#'
#' @param numiter The number of updates to perform.
#'
#' @param \dots Additional arguments passed to
#'   \code{\link{fit_poisson_nmf}}.
#'
#' @seealso \code{\link{fit_poisson_nmf}}
#' 
#' @examples
#' \dontrun{
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
#' dat <- simulate_multinom_gene_data(175,1200,k = 3)
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
#' }
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
    return(project_poisson_nmf(newdata,object$F,numiter,...)$L)
}

#' @rdname predict
#' 
#' @importFrom stats predict
#' 
#' @method predict multinom_topic_model_fit
#'
#' @export
#' 
predict.multinom_topic_model_fit <- function (object, newdata,
                                              numiter = 20, ...) {
  if (missing(newdata))
    return(object$L)
  else
    return(poisson2multinom(project_poisson_nmf(newdata,
                                                multinom2poisson(object)$F,
                                                numiter,...))$L)
}

# Fit a Poisson non-negative factorization to X in which the factors
# matrix is fixed to F.
project_poisson_nmf <- function (X, F, numiter, ...) {

  # Verify that fit0 and X are compatible.
  n <- nrow(X)
  k <- ncol(F)
  verify.count.matrix(X)
  if (nrow(F) != ncol(X))
    stop("fit$F and X do not match; the number of rows of fit$F should ",
         "equal the number of columns of X")

  # Fit a Poisson NMF to the data, X, with F fixed.
  L0  <- matrix(1/k,n,k)
  fit <- init_poisson_nmf(X,F = F,L = L0)
  fit <- fit_poisson_nmf(X,fit0 = fit,numiter = numiter,...)
  return(fit)
}
