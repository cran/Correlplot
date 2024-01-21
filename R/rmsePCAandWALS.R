rmsePCAandWALS <- function(R,output,digits=4,omit.diagonals=c(FALSE,FALSE,TRUE,TRUE)) {
  p <- ncol(R)
  W1 <- matrix(1,p,p)
  if(omit.diagonals[1]) rmse.pca <- rmse(R,output$Rhat.pca,per.variable = TRUE) else rmse.pca <- rmse(R,output$Rhat.pca,W1,per.variable = TRUE)
  if(omit.diagonals[2]) rmse.pca.adj <- rmse(R,output$Rhat.pca.adj,per.variable = TRUE) else rmse.pca.adj <- rmse(R,output$Rhat.pca.adj,W1,per.variable = TRUE)
  if(omit.diagonals[3]) rmse.wals <- rmse(R,output$Rhat.wals,per.variable = TRUE) else rmse.wals <- rmse(R,output$Rhat.wals,W1,per.variable = TRUE)
  if(omit.diagonals[4]) rmse.wals.adj <- rmse(R,output$Rhat.wals.adj,per.variable = TRUE) else rmse.wals.adj <- rmse(R,output$Rhat.wals.adj,W1,per.variable = TRUE)
  Results <- cbind(rmse.pca,rmse.pca.adj,rmse.wals,rmse.wals.adj)
  colnames(Results) <- c("PCA","PCA-A","WALS","WALS-A")
  if(omit.diagonals[1]) all.pca <- rmse(R,output$Rhat.pca) else all.pca <- rmse(R,output$Rhat.pca,W1)
  if(omit.diagonals[2]) all.pca.adj <- rmse(R,output$Rhat.pca.adj) else all.pca.adj <- rmse(R,output$Rhat.pca.adj,W1)
  if(omit.diagonals[3]) all.wals <- rmse(R,output$Rhat.wals) else all.wals <- rmse(R,output$Rhat.wals,W1)
  if(omit.diagonals[4]) all.wals.adj <- rmse(R,output$Rhat.wals.adj) else all.wals.adj <- rmse(R,output$Rhat.wals.adj,W1)
  allrow <- c(all.pca,all.pca.adj,all.wals,all.wals.adj)
  Results <- rbind(Results,allrow)
  Results <- round(Results,digits)
  rownames(Results)[nrow(Results)] <- "All"
  return(Results)
}
