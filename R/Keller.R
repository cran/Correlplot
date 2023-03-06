Keller <- function(R, eps = 1e-6, nd = 2, itmax = 10) {
  out.kel <- eigen(R)
  V <- out.kel$vectors
  Dl <- diag(out.kel$values)
  Rhat <- V[,1:nd]%*%Dl[1:nd,1:nd]%*%t(V[,1:nd])
  e <- diag(R) - diag(Rhat)
  esq <- sum(e*e)
  Radj <- R
  diag(Radj) <- diag(Rhat)
  i <- 1
  cat(i,esq,"\n")
  while (esq > eps & i < itmax) {
    out.kel <- eigen(Radj)
    V <- out.kel$vectors
    Dl <- diag(out.kel$values)
    Rhat <- V[,1:nd]%*%Dl[1:nd,1:nd]%*%t(V[,1:nd])
    e <- diag(Radj) - diag(Rhat)
    esq <- sum(e*e)
    diag(Radj) <- diag(Rhat)
    i <- i + 1
    cat(i,esq,"\n")
  }
  return(Rhat)
}
