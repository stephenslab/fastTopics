# Verify that x is a minimally valid non-negative matrix.
verify.matrix <- function (x, arg.name = deparse(substitute(x))) {
  arg.name <- sprintf("\"%s\"",arg.name)
  msg <- paste("Input argument",arg.name,"should be a non-negative,",
               "numeric matrix (a \"matrix\" or a \"dgCMatrix\") with at",
               "least two 2 rows and at least 2 columns, and all entries",
               "should be finite and non-missing")
  if (!((is.matrix(x) & is.numeric(x)) |
        inherits(x,"dgCMatrix")))
    stop(msg)
  else if (nrow(x) < 2 |
           ncol(x) < 2 |
           any(x < 0)  |
           any(is.infinite(x)) |
           any(is.na(x)))
    stop(msg)
  return(TRUE)
}

# Verify that x is a list with elements F and L.
verify.fit <- function (x, arg.name = deparse(substitute(x))) {
  arg.name <- sprintf("\"%s\"",arg.name)
  msg <- paste("Input argument",arg.name,"should be a list containing",
               "named elements \"F\" and \"L\"")
  if (!is.list(x))
    stop(msg)
  else if (!all(is.element(c("F","L"),names(x))))
    stop(msg)
  return(TRUE)
}

# Apply operation f to all nonzeros of a sparse matrix.
apply.nonzeros <- function (X, f) {
  d <- summary(X)
  return(sparseMatrix(i = d$i,j = d$j,x = f(d$x),dims = dim(X)))
}

# scale.cols(A,b) scales each column A[,i] by b[i].
#
# scale.cols(A) scales each column of A so that the entries of each
# column sum to 1; this is the same as scale.cols(A,1/colSums(A)).
scale.cols <- function (A, b) {
  if (missing(b))
    b <- 1/colSums(A)
  return(t(t(A) * b))
}

# Scale each row of A so that the entries of each row sum to 1.
normalize.rows <- function (A)
  A / rowSums(A)

# Scale each column of A so that the entries of each column sum to 1.
normalize.cols <- function (A)
  t(t(A) / colSums(A))

# Rescale the factors (F) and loadings (F) with the property that the
# matrix reconstruction L*F' remains the same after rescaling;
# specifically, rescale the columns of F and L so that, for each k,
# column k of F has the same mean as column k of L.
rescale.factors <- function (F, L) {
  d <- sqrt(colMeans(L)/colMeans(F))
  return(list(F = F/d,L = d*L))
}

# trcrossprod(A,B) returns trace(t(A) %*% B), where A and B are
# matrices with the same number of rows and column, and trace(X) =
# sum(diag(X)).
trcrossprod <- function (A, B)
  sum(A * B)
