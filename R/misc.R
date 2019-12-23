# Return TRUE if x is a finite scalar with no missing entries.
is.scalar <- function (x)
  is.numeric(x) &
  length(x) == 1 &
  all(!is.na(x)) &
  all(is.finite(x))

# Verify that x is non-negative matrix.
verify.nonnegative.matrix <- function (x, arg.name = deparse(substitute(x))) {
  arg.name <- sprintf("\"%s\"",arg.name)
  msg <- paste("Input argument",arg.name,"should be a non-negative,",
               "numeric matrix (a \"matrix\" or a \"dgCMatrix\"), and",
               "all entries should be finite and non-missing")
  if (!((is.matrix(x) & is.numeric(x)) | inherits(x,"dgCMatrix")))
    stop(msg)
  else if (any(x < 0) | any(is.infinite(x)) | any(is.na(x)))
    stop(msg)
  return(TRUE)
}

# Verify that x is a valid count matrix.
verify.count.matrix <- function (x, arg.name = deparse(substitute(x))) {
  verify.nonnegative.matrix(x,arg.name)
  arg.name <- sprintf("\"%s\"",arg.name)
  if (!(nrow(x) > 1 & ncol(x) > 1))
    stop(paste("Input matrix",arg.name,"should have at least 2 rows",
               "and 2 columns"))
  return(TRUE)
}

# Verify that x is a valid topic model fit or non-negative matrix
# factorization.
verify.fit <- function (x, arg.name = deparse(substitute(x))) {
  arg.name.F <- paste0(arg.name,"$F")
  arg.name.L <- paste0(arg.name,"$L")
  arg.name   <- sprintf("\"%s\"",arg.name)
  msg        <- paste("Input argument",arg.name,"should be a list containing",
                      "non-negative matrices \"F\" and \"L\"")
  if (!is.list(x))
    stop(msg)
  else if (!all(is.element(c("F","L"),names(x))))
    stop(msg)
  verify.nonnegative.matrix(x$F,arg.name.F)
  verify.nonnegative.matrix(x$L,arg.name.L)
  if (ncol(x$F) != ncol(x$L))
    stop(paste("Input matrices",arg.name.F,"and",arg.name.L,"should have",
               "the same number of columns"))
  return(TRUE)
}

# Verify that x is a valid count matrix and "fit" is a valid topic model
# fit or non-negative matrix factorization.
verify.fit.and.count.matrix <-
    function (x, fit,
              arg.name.x   = deparse(substitute(x)),
              arg.name.fit = deparse(substitute(fit))) {
  verify.count.matrix(x,arg.name.x)
  verify.fit(fit,arg.name.fit)
  arg.name.x <- sprintf("\"%s\"",arg.name.x)
  arg.name.F <- sprintf("\"%s$F\"",arg.name.fit)
  arg.name.L <- sprintf("\"%s$L\"",arg.name.fit)
  if (!(nrow(fit$L) == nrow(x) & nrow(fit$F) == ncol(x)))
    stop(paste("Dimensions of input matrices",arg.name.x,",",arg.name.F,
               "and",arg.name.L,"do not agree"))
  return(TRUE)
}

# Apply operation f to all nonzeros of a sparse matrix.
apply.nonzeros <- function (X, f) {
  d <- summary(X)
  return(sparseMatrix(i = d$i,j = d$j,x = f(d$x),dims = dim(X)))
}

# scale.cols(A,b) scales each column A[,i] by b[i].
scale.cols <- function (A, b)
  t(t(A) * b)

# Scale each row of A so that the entries of each row sum to 1.
normalize.rows <- function (A)
  A / rowSums(A)

# Scale each column of A so that the entries of each column sum to 1.
normalize.cols <- function (A)
  t(t(A) / colSums(A))

# Scale each row of A so that the large entry in each row is 1.
normalize.rows.by.max <- function (A)
  A / apply(A,1,max)

# Rescale the factors (F) and loadings (F) with the property that the
# matrix reconstruction L*F' remains the same after rescaling;
# specifically, rescale the columns of F and L so that, for each k,
# column k of F has the same mean as column k of L.
rescale.factors <- function (F, L) {
  d <- sqrt(colMeans(L)/colMeans(F))
  return(list(F = scale.cols(F,d),
              L = scale.cols(L,1/d)))
}

# For a Poisson non-negative matrix factorization with rank = 1, the
# maximum-likelihood estimate (MLE) has a closed-form (up to a scaling
# factor); this function returns the MLE subject to the constraint
# that mean(F) = mean(L).
fit_pnmf_rank1 <- function (X)
  list(F = matrix(colMeans(X)),L = matrix(rowMeans()))
