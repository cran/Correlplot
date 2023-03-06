wOne <-
function(x, w) {
  n <- nrow(x)
  m <- ncol(x)
  d <- sum(w * x) / sum(w)
  p <- rep(1, n)
  q <- rep(1, m)
  delta <- d * matrix(1, n, m)
  return(list(p = p, q = q, delta = delta))
}
