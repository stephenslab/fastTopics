# Verify that x is a vector with positive entries.
verify.positive.vector <- function (x, arg.name = deparse(substitute(x))) {
  arg.name <- sprintf("\"%s\"",arg.name)
  msg <- paste("Input argument",arg.name,"should be a numeric vector in",
               "which all entries are finite, non-missing and positive")
  if (!is.numeric(x))
    stop(msg)
  else if (any(x <= 0) | any(is.infinite(x)) | any(is.na(x)))
    stop(msg)
  return(TRUE)
}

# Verify that x is non-negative matrix.
verify.nonnegative.matrix <- function (x, arg.name = deparse(substitute(x))) {
  arg.name <- sprintf("\"%s\"",arg.name)
  msg <- paste("Input argument",arg.name,"should be a non-negative,",
               "numeric matrix (a \"matrix\" or a \"dgCMatrix\"), and",
               "all entries should be finite and non-missing")
  if (!((is.matrix(x) & is.numeric(x)) | is.sparse.matrix(x)))
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

# Verify that x is a valid multinomial topic model fit or non-negative
# matrix factorization.
verify.fit <- function (x, arg.name = deparse(substitute(x))) {
  arg.name.F <- paste0(arg.name,"$F")
  arg.name.L <- paste0(arg.name,"$L")
  arg.name.s <- paste0(arg.name,"$s")
  arg.name   <- sprintf("\"%s\"",arg.name)
  msg        <- paste("Input argument",arg.name,"should be a list containing",
                      "non-negative matrices \"F\" and \"L\"")
  if (!is.list(x))
    stop(msg)
  else if (!all(is.element(c("F","L"),names(x))))
    stop(msg)
  verify.nonnegative.matrix(x$F,arg.name.F)
  verify.nonnegative.matrix(x$L,arg.name.L)
  arg.name.F <- sprintf("\"%s\"",arg.name.F)
  arg.name.L <- sprintf("\"%s\"",arg.name.L)
  if (ncol(x$F) != ncol(x$L))
    stop(paste("Input matrices",arg.name.F,"and",arg.name.L,"should have",
               "the same number of columns"))
  if (is.element("s",names(x))) {
    
    # Check the vector of "scale factors", s.
    verify.positive.vector(x$s,arg.name.s)
    arg.name.s <- sprintf("\"%s\"",arg.name.s)
    if (length(x$s) != nrow(x$L))
      stop(paste("The length of input vector",arg.name.s,"should equal the",
           "number of rows in",arg.name.L))
  }
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

# Return TRUE if x is a finite scalar with no missing entries.
is.scalar <- function (x)
  is.numeric(x) &
  length(x) == 1 &
  all(!is.na(x)) &
  all(is.finite(x))
