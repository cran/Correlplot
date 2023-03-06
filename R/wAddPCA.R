wAddPCA <- function (x, w = matrix(1, nrow(x), ncol(x)), p = 2, add = "all", 
    bnd = "opt", itmaxout = 1000, itmaxin = 1000, epsout = 1e-06, 
    epsin = 1e-06, verboseout = TRUE, verbosein = FALSE) 
{
    if ((add == "row") && isSymmetric(x)) {
        warning("row option is problematic for symmetric data")
    }
    if ((add == "col") && isSymmetric(x)) {
        warning("col option is problematic for symmetric data")
    }
    delta <- 0
    itel <- 1
    s <- svd(x)
    if(p==1) {
      zold <- tcrossprod(s$u[,1]*s$d[1], s$v[,1])
    } else {
      zold <- tcrossprod(s$u[, 1:p] %*% diag(s$d[1:p]), s$v[, 1:p])  
    }
    fold <- sum(w * (x - zold)^2)
    repeat {
        h <- wPCA(x - delta, w, p = p, zold = zold, bnd = bnd, 
            itmax = itmaxin, eps = epsin, verbose = verbosein)
        fmid <- h$f
        znew <- h$z
        if (add == "nul") {
            g <- wNul(x - znew, w)
            delta <- g$delta
        }
        if (add == "one") {
            g <- wOne(x - znew, w)
            delta <- g$delta
        }
        if (add == "row") {
            g <- wRow(x - znew, w)
            delta <- g$delta
        }
        if (add == "col") {
            g <- wCol(x - znew, w)
            delta <- g$delta
        }
        if (add == "all") {
            g <- wAll(x - znew, w)
            delta <- g$delta
        }
        fnew <- sum(w * (x - (delta + znew))^2)
        if (verboseout) {
            cat("iteration ", formatC(itel, digits = 0, width = 3, 
                format = "d"), " fold", formatC(fold, digits = 6, 
                width = 10, format = "f"), " fmid", formatC(fmid, 
                digits = 6, width = 10, format = "f"), " fnew", 
                formatC(fnew, digits = 6, width = 10, format = "f"), 
                "\n")
        }
        if (((fold - fnew) < epsout) || (itel == itmaxout)) {
            if(itel == itmaxout) {
              cat("Maximum number of iterations (",itmaxout,") reached.\n")
              cat("The algorithm may not have converged.\n")
            }
            break
        }
        itel <- itel + 1
        fold <- fnew
        zold <- znew
    }
    return(list(a = h$a, b = h$b, z = h$z, f = h$f, u = h$u, 
        v = h$v, p = g$p, q = g$q, itel = itel, delta = delta, 
        y = h$z + delta))
}