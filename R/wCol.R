wCol <-
function(x, w) {
  n <- nrow(x)
  p <- rep(0, n)
  q <- colSums(x * w) / colSums(w)
  delta <- outer(p, q, "+")
  return(list(p = p, q = q, delta = delta))
}
