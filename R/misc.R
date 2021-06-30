# Compute two-tailed p-value from z-score.
#
#' @importFrom stats pnorm
pfromz <- function (z)
  2*pnorm(-abs(z))

# Compute log10 two-tailed p-value from z-score.
#
#' @importFrom stats pnorm
lpfromz <- function (z)
  (log(2) + pnorm(-abs(z),log.p = TRUE))/log(10)

# Return true if x is a compressed, sparse, column-oriented numeric
# matrix.
is.sparse.matrix <- function (x)
  inherits(x,"dgCMatrix")

# Efficiently extract the nonzero elements from column j of sparse
# matrix A (a member of class "dgCMatrix"). Output "x" contains the
# nonzero values, and output "i" contains the
get.nonzeros <- function (A, j)
  list(x = A[,j,drop = FALSE]@x,i = A[,j,drop = FALSE]@i + 1)

# Check if the matrix contains one or more all-zero columns.
#
#' @importFrom Matrix colSums
any.allzero.cols <- function (X)
  any(colSums(X > 0) == 0)

# Filter out all-zero columns from the matrix.
#
#' @importFrom Matrix colSums
remove.allzero.cols <- function (X)
  X[,colSums(X > 0) >= 1]

# Apply operation f to all nonzeros of a sparse matrix.
#
#' @importFrom Matrix sparseMatrix
#' 
apply.nonzeros <- function (X, f) {
  d <- summary(X)
  return(sparseMatrix(i = d$i,j = d$j,x = f(d$x),dims = dim(X)))
}

# Compute X/(crossprod(A,B) + e) efficiently when X is a sparse
# matrix.
#
#' @importFrom Matrix sparseMatrix
#' @importFrom Rcpp evalCpp
#' 
x_over_tcrossprod <- function (X, A, B, e) {
  d <- summary(X)
  y <- drop(x_over_crossprod_rcpp(d$i - 1,d$j - 1,d$x,A,B,e))
  return(sparseMatrix(i = d$i,j = d$j,x = y,dims = dim(X)))
}

# Return an m x n matrix rbind(x,...,x), in which length(x) = m.
repmat <- function (x, n)
  matrix(x,n,length(x),byrow = TRUE)

# scale.cols(A,b) scales each column A[,i] by b[i].
scale.cols <- function (A, b)
  t(t(A) * b)

# Scale each row of A so that the entries of each row sum to 1.
#
#' @importFrom Matrix rowSums
#'
normalize.rows <- function (A)
  A / rowSums(A)

# Scale each column of A so that the entries of each column sum to 1.
#
#' @importFrom Matrix colSums
#'
normalize.cols <- function (A)
  t(t(A) / colSums(A))

# Scale each row of A so that the large entry in each row is 1.
normalize.rows.by.max <- function (A) {
  if (!is.matrix(A))
    stop("Input argument \"A\" should be a matrix")
  return(A / apply(A,1,max))
}
  
# Rescale the factors (F) and loadings (F) with the property that the
# matrix reconstruction L*F' remains the same after rescaling;
# specifically, rescale the columns of F and L so that, for each k,
# column k of F has the same mean as column k of L.
#
#' @importFrom Matrix colMeans
#'
rescale.factors <- function (F, L) {
  d <- sqrt(colMeans(L)/colMeans(F))
  return(list(F = scale.cols(F,d),
              L = scale.cols(L,1/d)))
}

# This does the same thing as the "rand" function in MATLAB.
#
#' @importFrom stats runif
#'
rand <- function (n, m, min = 0, max = 1) 
  matrix(runif(n*m,min,max),n,m)

# Initialize RcppParallel multithreading using a pre-specified number
# of threads, or using the default number of threads when "n" is NA.
#
#' @importFrom RcppParallel setThreadOptions
#' @importFrom RcppParallel defaultNumThreads
#'
initialize.multithreading <- function (n) {
  if (is.na(n)) {
    setThreadOptions()
    n <- defaultNumThreads()
  } else
    setThreadOptions(numThreads = n)
  if (n > 1)
    message(sprintf("Using %d RcppParallel threads.",n))
  return(n)
}

# For a Poisson non-negative matrix factorization with rank = 1, the
# maximum-likelihood estimate (MLE) has a closed-form solution (up to
# a scaling factor); this function returns the MLE subject to the
# constraint that mean(F) = mean(L).
#
#' @importFrom Matrix rowMeans
#' @importFrom Matrix colMeans
#'
fit_pnmf_rank1 <- function (X)
  list(F = matrix(colMeans(X)),
       L = matrix(rowMeans(X)))

# Compute the highest posterior density (HPD) interval from a vector
# of random draws from the distribution. See Chen & Shao (1999) for
# background on HPD intervals.
hpd <- function (x, conf.level = 0.95) {
  n <- length(x)
  m <- round(n*(1 - conf.level))
  x <- sort(x)
  y <- x[seq(n-m+1,n)] - x[seq(1,m)]
  i <- which.min(y)
  return(c(x[i],x[n-m+i]))
}
