FitRwithPCAandWALS <- function(R,nd=2,itmaxout=10000,
                               itmaxin=10000,eps=1e-8) {
  #
  #
  #
  p <- ncol(R)
  W1 <- matrix(1,ncol=p,nrow=p)
  W0 <- matrix(1,ncol=p,nrow=p)
  diag(W0) <- 0
  #
  # PCA
  #
  out.sd <- eig(R)
  V  <- out.sd$V[,1:nd]
  Ds <- sqrt(out.sd$D[1:nd,1:nd])
  Gp <- V%*%Ds
  Rhat.pca <- Gp%*%t(Gp)
  #
  # PCA with adjustment: WALS with diagonal turned on; 
  # and delta adjustment
  #
  Results <- wAddPCA(R, W1, add = "one", verboseout = FALSE, 
                   epsout = eps, itmaxout = itmaxout, 
                   itmaxin = itmaxin, p = nd)
  Rhat.pca.adj <- Results$a%*%t(Results$b) + Results$delta
  #
  # WALS with diagonal turned OFF and no delta adjustment
  #
  Results <- wAddPCA(R, W0, add = "nul", verboseout = FALSE, 
                   epsout = eps, itmaxout = itmaxout, 
                   itmaxin = itmaxin, p = nd)
  Rhat.wals <- Results$a%*%t(Results$b)
  #
  # WALS with delta adjustment
  #
  Results <- wAddPCA(R, W0, add = "one", verboseout = FALSE, 
                   epsout = eps, itmaxout = itmaxout, 
                   itmaxin = itmaxin, p = nd)
  Rhat.wals.adj <- Results$a%*%t(Results$b) + Results$delta
  return(list(Rhat.pca=Rhat.pca,Rhat.pca.adj=Rhat.pca.adj,
              Rhat.wals=Rhat.wals,Rhat.wals.adj=Rhat.wals.adj))
}
