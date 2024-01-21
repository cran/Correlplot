eig <- function(X) {
  out <- eigen(X)
  V <- out$vectors
  D <- diag(out$values)
  return(list(V=V,D=D))
}
