# Verify that x is a minimally valid non-negative matrix.
verify.matrix <- function (x, arg.name = deparse(substitute(x))) {
  arg.name <- sprintf("\"%s\"",arg.name)
  msg <- paste("Input argument",arg.name,"should be a non-negative,",
               "numeric matrix with at least two 2 rows and at least",
               "2 columns, and all entries should be finite and non-missing")
  if (!((is.matrix(x) & is.numeric(x)) |
        inherits(x,"dMatrix")))
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

# scale.cols(A,b) scales each column A[,i] by b[i].
#
# scale.cols(A) scales each column of A so that the entries of each
# column sum to 1; this is the same as scale.cols(A,1/colSums(A)).
scale.cols <- function (A, b) {
  if (missing(b))
    b <- 1/colSums(A)
  return(t(t(A) * b))
}

