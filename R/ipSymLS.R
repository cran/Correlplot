ipSymLS <- function(target, w=matrix(1,dim(target)[1], dim(target)[2]),
                    ndim=2, init=FALSE, itmax=100, eps=1e-6, verbose=FALSE) {
  if(is.matrix(init)) 
    x <- init 
  else {
    z <- eigen(target); x <- z$vectors[,1:ndim]; x <- x%*%diag(sqrt(z$values[1:ndim]))
  }
  n <- dim(target)[1]; xx <- tcrossprod(x); oloss <- sum(w*(target-xx)^2); itel <- 1
  repeat{ 
    for(i in 1:n) {
      ai <- crossprod(x,w[i,]*x)-w[i,i]*outer(x[i,],x[i,])
      bi <- colSums((w[i,]*target[i,])*x)-w[i,i]*target[i,i]*x[i,]
      if(w[i,i]==0)
        x[i,] <- solve(ai,bi)
      else {
        li <- sum(x[i,]^2); zi <- x[i,]/sqrt(li); wi <- w[i,i]; ci <- target[i,i]
        a4 <- wi; a3 <- 0; a2 <- 2*(sum(zi*(ai%*%zi)) - (wi*ci)); a1 <- -4*sum(zi*bi)
        b1 <- a1/(2*a4); b2 <- a2/(3*wi); bb <- as.complex((b1^4)+2*(b1^2)*(b2^3))
        s <- (b1^2)+(b2^3)+sqrt(bb); kf <- s^(1/3) + (b2^2)*s^(-1/3)+b2; li <- Re(-b1/kf)
        x[i,] <- li*zi; aa <- 2*(li^2)*ai; bb <- 2*li*bi; ei <- eigen(aa)
        l <- ei$values; k <- ei$vector; kb <- drop(crossprod(k,bb)); ml <- min(l)
        ul <- ml-sqrt(sum(kb^2)); uu <- ml-sqrt(sum(kb[which(l==ml)]^2)); u <- uu
        fu <- function(u)sum((kb/(l-u))^2)-1
        gu <- function(u)2*sum((kb^2)/((l-u)^3))
        repeat{
          if(abs(fu(u)) < 1e-6) break()
          u <- u-fu(u)/gu(u)
        }
        x[i,] <- li*(k%*%(kb/(l-u)))
      }
   }
  xx <- tcrossprod(x); nloss <- sum(w*(target-xx)^2)
  if(verbose) cat("Iteration:",formatC(itel,digits=6,width=6),
                  "Previous Loss:",formatC(oloss,digits=6,width=12,
                                           format="f"),
                  "Current Loss:",formatC(nloss,digits=6,width=12,format="f"),
                  "\n")
  if(((oloss-nloss) < eps)) break()
  if((itel==itmax)) {
      cat(paste("Maximum number of inner iterations (",itmax,") reached.\n",sep=""))
      cat("Algorithm ipSymLS may not have converged.\n") 
     break()
  }
  oloss <- nloss;
  itel <- itel + 1
}
return(x)
}
