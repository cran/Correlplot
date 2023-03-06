rmsePCAandWALS <- function(R,output,digits=4,diagonals=c(FALSE,FALSE,TRUE,TRUE)) {
  rmse.pca       <- rmse(R,output$Rhat.pca,omit.diagonal=diagonals[1],
                         per.variable = TRUE)
  rmse.pca.adj   <- rmse(R,output$Rhat.pca.adj,omit.diagonal=diagonals[2],
                         per.variable = TRUE)
  rmse.wals      <- rmse(R,output$Rhat.wals,omit.diagonal=diagonals[3],
                         per.variable = TRUE)
  rmse.wals.adj  <- rmse(R,output$Rhat.wals.adj,omit.diagonal=diagonals[4],
                         per.variable = TRUE)
  Results <- cbind(rmse.pca,rmse.pca.adj,rmse.wals,rmse.wals.adj)
  colnames(Results) <- c("PCA","PCA-A","WALS","WALS-A")  
  all.pca      <- rmse(R,output$Rhat.pca,omit.diagonal=diagonals[1])
  all.pca.adj  <- rmse(R,output$Rhat.pca.adj,omit.diagonal=diagonals[2])    
  all.wals     <- rmse(R,output$Rhat.wals,omit.diagonal=diagonals[3])
  all.wals.adj <- rmse(R,output$Rhat.wals.adj,omit.diagonal=diagonals[4])
  allrow <- c(all.pca,all.pca.adj,all.wals,all.wals.adj)
  Results <- rbind(Results,allrow)
  Results <- round(Results,digits)
  rownames(Results)[nrow(Results)] <- "All"
  return(Results)
}
