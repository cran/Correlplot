FitRDeltaQSym <- function(R,W=NULL,nd=2,eps=1e-10,delta=0,
                          q=colMeans(R),itmax.inner=1000,
                          itmax.outer=1000,verbose=FALSE) {
  p <- ncol(R)
  W <- W/sum(W)
  J <- matrix(1,nrow=p,ncol=p)
  if(is.null(W)) {
    W <- matrix(1,p,p)
    diag(W) <- 0
    W <- W/sum(W)
  }
  outer.ssq <- 10
  itel   <- 1
  while(outer.ssq > eps & itel < itmax.outer) {
    inner.ssq <- 10
    jtel <- 1
    while(inner.ssq > eps & jtel < itmax.inner) {
      Ra <- R - delta*J - 0.5*rep(1,p)%o%q - 0.5*q%o%rep(1,p)
      C <- ipSymLS(Ra,W,ndim=nd,eps=eps,itmax=itmax.inner)
      Rh <- C%*%t(C)
      newdelta <- tr(R%*%W) - rep(1,p)%*%W%*%q - tr(Rh%*%W)
      newdelta <- as.vector(newdelta)
      aa <- (newdelta-delta)
      inner.ssq <- aa*aa
      if(verbose) cat(itel, jtel, inner.ssq,"\n")
      delta <- newdelta
      jtel <- jtel+1
    }
    if(jtel >= itmax.inner) {
      cat(paste("Maximum number of inner iterations (",itmax.inner,") reached.\n"))
      cat("The algorithm may not have converged.\n") 
    }
    Raa <- R - delta*J - Rh
    Dw.inv <- diag(1/rowSums(W))
    newq <- as.vector(Dw.inv%*%rowSums(W*Raa))
    change.q <- newq-q
    outer.ssq <- sum(change.q^2)
    if(verbose) cat(itel,outer.ssq,"\n")
    q <- newq
    itel <- itel + 1
  }
  if(itel >= itmax.outer) {
    cat(paste("Maximum number of outer iterations (",itmax.inner,") reached.\n"))
    cat("The algorithm may not have converged.\n") 
  }
  if(verbose) {
    cat("final delta:",delta,"\n")
    cat("final q:\n") 
    print(q)
  }
  Rhat <- delta*J + rep(1,p)%o%q + C%*%t(C)
  rmse.a <- rmse(R,Rhat)
  return(list(delta=delta,Rhat=Rhat,C=C,rmse=rmse.a,q=q))
}
