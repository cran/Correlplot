rmse <- function(R,Rhat,option="offd",verbose=FALSE) {
  #
  # calculate root mean squared error of approximation
  #
  p <- ncol(R)
  if(verbose) cat(p,"variables\n")
  E <- R - Rhat
  if(option=="full") {
    sse <- sum(E^2)
    mse <- sse/(p*p)
    y <- sqrt(mse)
    if(verbose) cat("rmse (full matrix) = ",y,"\n")
  } else if (option=="offd") {
    El <- E[lower.tri(E)]
    sse <- sum(El*El)
    mse <- sse/(0.5*p*(p-1))
    y <- sqrt(mse)
    if(verbose) cat("rmse (off-diagonal) = ",y,"\n")
  } else {
    stop("rmse: unknown value for option.")
  }
  return(y)
}
