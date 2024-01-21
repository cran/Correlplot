rmse <- function(R, Rhat, 
                     W = matrix(1,nrow(R),ncol(R)) - diag(nrow(R)),
                     verbose = FALSE, per.variable = FALSE)  { 
  p <- ncol(R)
  rmse.per.var <- rep(NA, p)
  names(rmse.per.var) <- colnames(R)
  rmse.overall <- NA
  E <- R - Rhat
  WE2 <- W * (E * E)
  num <- sum(WE2)
  den <- sum(W)
  rmse.overall <- sqrt(num/den)
  if(per.variable) {
    for(i in 1:p) {
      numi <- sum(WE2[i,]) + sum(WE2[,i]) - WE2[i,i]
      deni <- sum(W[i,]) + sum(W[,i]) - W[i,i]
      rmse.per.var[i] <- sqrt(numi/deni)
    }
  }
  if (verbose) {
    cat(p, "variables\n")
    cat("weight matrix:\n")
    print(W)
    cat("RMSE = ",rmse.overall,"\n")
    if (per.variable) {
      print(rmse.per.var)
    }
  }
  if(!per.variable) {
    return(rmse.overall)
  } else {
    return(rmse.per.var)
  }
}
