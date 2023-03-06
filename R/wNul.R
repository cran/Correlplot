wNul <-
function(x, w) {
  n <- nrow(x)
  m <- ncol(x)
  p <- rep(0, n)
  q <- rep(0, m)
  delta <- outer(p, q, "+")
  return(list(p = p, q = q, delta = delta))
}
