correlogram <- function (R, labs = colnames(R), ifun = "cos", cex = 1, 
                            main = "", ntrials = 50, 
                            xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2), 
                            pos = NULL, ...) 
{
  theta <- fit_angles(R, ifun = ifun, ntrials = ntrials)
  X <- matrix(c(cos(theta), sin(theta)), ncol = 2)
  if (is.null(labs)) 
    labs <- 1:ncol(R)
  opar <- par(pty = "s", bty = "n", xaxt = "n", 
              yaxt = "n")
  plot(0, 0, type = "n", xlab = "", ylab = "", 
       xlim = xlim, ylim = ylim, main = main, ...)
  #  points(X[, 1], X[, 2], pch = 19)
  if (is.null(pos)) textxy(X[, 1], X[, 2], labs, cex = cex) else
    text(X[,1],X[,2],labs, pos = pos, cex = cex)
  #  arrows(0, 0, X[, 1], X[, 2], length = 0)
  arrows(0, 0, X[, 1], X[, 2], length = 0.1, lwd = 1, 
         angle = 10, col="blue")
  circle()
  par(opar)
  return(theta)
}
