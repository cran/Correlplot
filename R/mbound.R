mbound <-
function (w) {
  n <- nrow (w)
  m <- ncol (w)
  g1 <- matrix(0, n * m, n)
  g2 <- matrix(0, n * m, m)
  k <- 1
  for (j in 1:m) {
    for (i in 1:n) {
      if (w[i, j] == 0) {
        next
      }
      g1[k, i] <- 1
      g2[k, j] <- 1
      k <- k + 1
    }
  }
  g <- cbind (g1, g2)[1:(k - 1), ]
  ww <- as.vector(w)
  f <- log (ww[ww > 0])
  h <- lsi (g, f, g, f)
  u <- exp (h[1:n])
  v <- exp (h[n + (1:m)])
  return (list (u = u, v = v))
}
