wRow <-
function(x, w) {
  m <- ncol(x)
  p <- rowSums(x * w) / rowSums(w)
  q <- rep(0, m)
  delta <- outer(p, q, "+")
  return(list(p = p, q = q, delta = delta))
}
