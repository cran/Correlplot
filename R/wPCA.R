wPCA <- function (x, w, p = 2, zold = lfunc(x, p), bnd = "opt", itmax = 1000, 
    eps = 1e-06, verbose = TRUE) 
{
    n <- nrow(x)
    m <- ncol(x)
    if (bnd == "all") {
        u <- rep(1, n)
        v <- rep(max(w), m)
    }
    if (bnd == "row") {
        u <- apply(w, 1, max)
        v <- rep(1, m)
    }
    if (bnd == "col") {
        u <- rep(1, n)
        v <- apply(w, 2, max)
    }
    if (bnd == "opt") {
        uv <- mbound(w)
        u <- uv$u
        v <- uv$v
    }
    fold <- sum(w * (x - zold)^2)
    itel <- 1
    repeat {
        z1 <- gfunc(zold, x, u, v, w)
        z2 <- lfunc(z1, p)
        znew <- sfunc(z2, u, v)
        fnew <- sum(w * (x - znew)^2)
        if (verbose) {
            cat("iteration ", formatC(itel, digits = 0, width = 3, 
                format = "d"), " fold", formatC(fold, digits = 6, 
                width = 10, format = "f"), " fnew", formatC(fnew, 
                digits = 6, width = 10, format = "f"), "\n")
        }
        if (((fold - fnew) < eps) || (itel == itmax)) {
          if(itel==itmax) {
            cat("Maximum number of inner iterations reached.\n")
          }
          break
        } 
        itel <- itel + 1
        zold <- znew
        fold <- fnew
    }
    snew <- svd(znew)
    if(p==1) {
      b <- snew$v[,1]*snew$d[1]
    } else {
      b <- snew$v[, 1:p] %*% diag(snew$d[1:p])
    }
    return(list(a = snew$u[, 1:p], b = b, 
        z = znew, f = fnew, itel = itel, u = u, v = v))
}