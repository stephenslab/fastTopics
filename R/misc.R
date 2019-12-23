# Apply operation f to all nonzeros of a sparse matrix.
#
#' @importFrom Matrix sparseMatrix
#' 
apply.nonzeros <- function (X, f) {
  d <- summary(X)
  return(sparseMatrix(i = d$i,j = d$j,x = f(d$x),dims = dim(X)))
}

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

# For a Poisson non-negative matrix factorization with rank = 1, the
# maximum-likelihood estimate (MLE) has a closed-form (up to a scaling
# factor); this function returns the MLE subject to the constraint
# that mean(F) = mean(L).
#
#' @importFrom Matrix rowMeans
#' @importFrom Matrix colMeans
#'
fit_pnmf_rank1 <- function (X)
  list(F = matrix(colMeans(X)),L = matrix(rowMeans(X)))
