# Much of the code contained here is based on the "TopicScore" package
# source code developed by Minzhe Wang and Tracy Ke, distributed under
# the MIT license.

# Estimates the word-topic matrix (A) from the word-document matrix
# (X) using the Topic SCORE algorithm.
#
# The inputs are: X, the n x m counts matrix (may be sparse or dense);
# k, the number of topics; k0, the number of greedy search steps to
# use in Vertex Hunting; m, the number of centers in the k-means step
# of Vertex Hunting; and Mquantile, the percentage of the quantile of
# the diagonal entries of matrix M, which is used to upper truncate
# the diagonal entries of matrix M.  When it is zero, it will
# degenerate the case when there is no normalization. When it's 1,
# there is no truncation.
#
# The return value is an m x k word-topic matrix.
#
#' @importFrom irlba irlba
#' 
topic_score <- function (X, k, k0 = ceiling(1.5*k), m = 3*k, Mquantile = 0) {

  # M0 = D0*1/n, where D0(j,i) is the expected frequency of word j
  # in document i.
  M0 <- colMeans(X)
  M0 <- pmin(M0,quantile(M0,Mquantile))
  
  # Compute the k right singular vectors of the normalized counts matrix.
  X <- scale.cols(X,1/sqrt(M0))
  V <- irlba(X,k)$v
  
  # Step 1: Recover the left-scaling matrix (LSM).
  v1 <- abs(V[,1])
  R  <- V[,-1]/v1
  
  # Step 2: Perform "Vertex Hunting".
  V <- vertex_hunting(R,k0,m)
  
  # Step 3: Recover the normalized topic matrix (NTM).
  P <- cbind(R,1) %*% solve(cbind(V,1))
  P <- pmax(P,0)
  P <- P / rowSums(P)
  
  # Step 4: Recover the unscaled topic matrix.
  A <- sqrt(M0)*v1*P

  # Step 5: Return the scaled topic matrix.
  return(normalize.cols(A))
}

# The Vertex Hunting algorithm for Topic-SCORE.  It finds a simplex
# with k vertices that best approximates the given p data points in a
# (k-1) dimensional space.
#
# The inputs are: R, the p x k-1 data matrix, with each row being a
# data point; k0, the number of greedy search steps; and m, the number of
# centers in the k-means step.
#
# The output is the k x k-1 vertices matrix, with each row being a
# vertex in the found simplex.
#
#' @importFrom utils combn
#' @importFrom stats kmeans
#' 
vertex_hunting <- function (R, k0, m) {
  k <- ncol(R) + 1
  
  # Step 2a.
  # X <- KMeans_rcpp(R,m,initializer = "kmeans++")$centroids
  X <- kmeans(R,m,iter.max = 100)$centers

  # Step 2b'.
  Y  <- tcrossprod(X)
  D  <- matrix(diag(Y),m,m)
  D  <- D + t(D) - 2*Y
  i  <- drop(arrayInd(which.max(D),dim(D)))
  X0 <- X[i,,drop = FALSE]
  X  <- X[-i,,drop = FALSE]
  if (k0 > 2) {
    for (j in 3:k0) {
      D  <- matrix(diag(tcrossprod(X)),j-1,nrow(X),byrow = TRUE)
      D  <- D - 2*tcrossprod(X0,X)
      i  <- which.max(colMeans(D))
      X0 <- rbind(X0,X[i,])
      X  <- X[-i,,drop = FALSE]
    }
    X <- X0
  }
  
  # Step 2b.
  B <- combn(1:k0,k)
  n <- ncol(B)
  v <- rep(0,n)
  for (i in 1:n)
    for (j in 1:k0)
      v[i] <- max(simplex_dist(X[j,],X[B[,i],,drop = FALSE]),v[i])
  i <- which.min(v)
  return(X[B[,i],])
}

# This function computes the shortest (Euclidean) distance between the
# given point (x) and any point in the simplex (V).
#
#' @importFrom quadprog solve.QP
#' 
simplex_dist <- function (x, V) {
  n  <- nrow(V)
  v  <- V[n,]
  A  <- cbind(diag(n-1),-1)
  VV <- A %*% V
  M  <- tcrossprod(VV)
  d  <- VV %*% (x - v)
  b0 <- rep(0,n)
  b0[n] <- -1
  f  <- solve.QP(M,d,A,b0)$value
  return(sqrt(max(sum((x - v)^2) + 2*f,0)))
}
