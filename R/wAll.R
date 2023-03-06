wAll <-
function(x, w) {
  n <- nrow(x)
  m <- ncol(x)
  b <- c(rowSums(w * x), colSums(w * x))
  a <- matrix(0, n + m, n + m)
  a[1:n, 1:n] <- diag(rowSums(w))
  a[n + 1:m, n + 1:m] <- diag(colSums(w))
  a[1:n, n + 1:m] <- w
  a[n + 1:m, 1:n] <- t(w)
  pq <- drop(ginv(a) %*% b)
  p <- pq[1:n]
  q <- pq[n + 1:m]
  delta <- outer(p, q, "+")
  return(list(p = p, q = q, delta = delta))
}
