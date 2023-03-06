rmse <- function(R,Rhat,W=matrix(1,nrow=nrow(R),ncol=ncol(R)),omit.diagonal=TRUE,verbose=FALSE,per.variable=FALSE) {
  #
  # calculate root mean squared error of approximation
  #
  p <- ncol(R)
  if(verbose) cat(p,"variables\n")
  if(omit.diagonal) {
     diag(W) <- 0
  }
  E  <- R - Rhat
  E2 <- W*(E*E)
  nonzero <- rowSums(W!=0)
  rmse.per.var <- sqrt(rowSums(E2)/nonzero)
  if(!omit.diagonal) {
    nelem <- sum(W!=0)
    E2 <- W*(E*E)
    sse <- sum(E2)
    mse <- sse/(nelem)
    y <- sqrt(mse)
    if(verbose) {
      cat("rmse (full matrix) = ",y,"\n")
      print(rmse.per.var)
    }
  } else {
    El <- E[lower.tri(E)]
    Wl <- W[lower.tri(W)]
    nelem <- sum(Wl!=0)
    El2 <- Wl*El*El
    sse <- sum(El2)
    mse <- sse/(nelem)
    y <- sqrt(mse)
    if(verbose) {
      cat("rmse (below-diagonal) = ",y,"\n")
      print(rmse.per.var)
    }
  }
  if(!per.variable) {
    return(y)
  } else {
    return(rmse.per.var)
  }
}
