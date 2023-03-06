lfunc <-
function (x, p) {
  l <- as.matrix(eigen (crossprod(x))$vectors[, 1:p])
  return (x %*% tcrossprod (l))
}
