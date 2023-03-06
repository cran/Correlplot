tr <-
function(X) 
{
  if (is.matrix(X)) {
    y <- sum(diag(X))
    n <- nrow(X)
    p <- ncol(X)
    if(n != p) {
      cat("warning: X is not square\n")
    }
  } else if (is.vector(X)) {
    cat("warning: X is not a matrix\n")
    y <- sum(X)
  }  else {
    cat("warning: X is nor matrix nor vector\n")
    y <- NA
  }
  return(y)
}
