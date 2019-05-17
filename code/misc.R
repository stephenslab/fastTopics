# Return the dot product of vectors x and y.
dot <- function (x, y)
  sum(x*y)

# Compute the value of the cost function for non-negative matrix
# factorization, in which matrix X is approximated by matrix AB = A *
# B. This is equivalent to the negative Poisson log-likelihood (after
# removing terms that do not depend on A and B).
cost <- function (X, AB, e = 1e-15)
  sum(AB - X*log(AB + e))

