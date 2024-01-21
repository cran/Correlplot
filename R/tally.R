tally <- function(G,adj=0,values=seq(-1,1,by=0.2),
                       pch=19,dotcolor="black",
                       cex=0.5,color.negative="red",
                       color.positive="blue") {
  if(is.vector(G)) {
    G <- matrix(G,ncol=2)
    Dg <- matrix(G%*%t(G),1,1)
  } else {
    G <- G[,1:2]
    Dg <- diag(diag(G%*%t(G)))
  }
  for(i in 1:nrow(G)) { 
    for(j in 1:length(values)) {
      DP <- (values[j]+adj)*G[i,]/Dg[i,i] # Dot positions
      points(DP[1],DP[2],pch=pch,cex=cex,col=dotcolor)
      if(j > 1) {
        DPold <- (values[j-1]+adj)*G[i,]/Dg[i,i]
        if((values[j-1]+adj) >= 0) {
          segments(DPold[1],DPold[2],DP[1],DP[2],lty="solid",col=color.positive)    
        } else {
          segments(DPold[1],DPold[2],DP[1],DP[2],lty="solid",col=color.negative)    
        }
      }
    }
  }
}
