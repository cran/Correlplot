## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
knitr::opts_chunk$set(fig.align = "center")
knitr::opts_chunk$set(fig.width = 6, fig.height = 6) 

## ----preinstall---------------------------------------------------------------
library(calibrate)
library(corrplot)
library(Correlplot)

## -----------------------------------------------------------------------------
data("Kernels")
X <- Kernels[Kernels$variety==1,]
X <- X[,-8]
head(X)

## -----------------------------------------------------------------------------
R <- cor(X)
round(R,digits=3)

## ----corrgram, fig.cap = "A corrgram of the wheat kernel data."---------------
corrplot(R, method="circle",type="lower")

## ----names, echo = FALSE------------------------------------------------------
vnames <- substr(colnames(R),1,4)
colnames(R) <- rownames(R) <- vnames
colnames(X) <- vnames

## ----correlogram, fig.cap = "The correlogram of the wheat kernel data."-------
theta.cos <- correlogram(R,xlim=c(-1.1,1.1),ylim=c(-1.1,1.1),main="CRG")

## -----------------------------------------------------------------------------
Rhat.cor <- angleToR(theta.cos)

## -----------------------------------------------------------------------------
rmse.crg <- rmse(R,Rhat.cor,omit.diagonal=FALSE)
rmse.crg

## ----linearcorrelogram, fig.cap="Linear correlogram of the wheat kernel data."----
theta.lin <- correlogram(R,ifun="lincos",labs=colnames(R),xlim=c(-1.1,1.1),
                         ylim=c(-1.1,1.1),main="CRG")

## -----------------------------------------------------------------------------
Rhat.corlin <- angleToR(theta.lin,ifun="lincos")
rmse.lin <- rmse(R,Rhat.corlin,omit.diagonal=FALSE)
rmse.lin

## ----pcabiplot, fig.cap = "PCA biplot of the wheat kernel data."--------------
n <- nrow(X)
Xt <- scale(X)/sqrt(n-1)
res.svd <- svd(Xt)
Fs <- sqrt(n)*res.svd$u # standardized principal components
Gp <- res.svd$v%*%diag(res.svd$d) # biplot coordinates for variables
bplot(Fs,Gp,colch=NA,collab=colnames(X),xlab="PC1",ylab="PC2",main="PCA")

## ----pcacorrelationcircle, fig.cap = "PCA biplot of the correlation matrix."----
bplot(Gp,Gp,colch=NA,rowch=NA,collab=colnames(X),xl=c(-1,1),
      yl=c(-1,1),main="PCA")
circle()

## -----------------------------------------------------------------------------
Rhat.pca <- Gp[,1:2]%*%t(Gp[,1:2])

## -----------------------------------------------------------------------------
rmse.pca <- rmse(R,Rhat.pca,omit.diagonal=FALSE)
rmse.pca

## -----------------------------------------------------------------------------
rmse(R,Rhat.pca,omit.diagonal=FALSE,per.variable=TRUE)

## ----mdsplot, fig.cap = "MDS map of the correlation matrix of the wheat kernel data."----
Di <- sqrt(2*(1-R))
out.mds <- cmdscale(Di,eig = TRUE)
Fp <- out.mds$points
opar <- par(bty = "l")
plot(Fp[,1],Fp[,2],asp=1,xlab="First principal axis",
     ylab="Second principal axis",main="MDS")
textxy(Fp[,1],Fp[,2],colnames(R),cex=0.75)
par(opar)

ii <- which(R < 0,arr.ind = TRUE)

for(i in 1:nrow(ii)) {
  segments(Fp[ii[i,1],1],Fp[ii[i,1],2],
           Fp[ii[i,2],1],Fp[ii[i,2],2],col="red",lty="dashed")
}

## -----------------------------------------------------------------------------
Dest <- as.matrix(dist(Fp[,1:2]))
Rhat.mds <- 1-0.5*Dest*Dest

## -----------------------------------------------------------------------------
rmse.mds <- rmse(R,Rhat.mds,omit.diagonal=FALSE)
rmse.mds

## ----pfa----------------------------------------------------------------------
out.pfa <- Correlplot::pfa(X,verbose=FALSE)
L <- out.pfa$La

## -----------------------------------------------------------------------------
Rhat.pfa <- L[,1:2]%*%t(L[,1:2])
rmse.pfa <- rmse(R,Rhat.pfa)
rmse.pfa

## ----pfabiplot, fig.cap = "PFA biplot of the correlation matrix of the wheat kernel data."----
bplot(L,L,pch=NA,xl=c(-1,1),yl=c(-1,1),
     xlab="Factor 1",ylab="Factor 2",main="PFA",rowch=NA,
     colch=NA)
circle()

## -----------------------------------------------------------------------------
diag(out.pfa$Psi)

## -----------------------------------------------------------------------------
W <- matrix(1,nrow(R),nrow(R))
diag(W) <- 0
Fp.als <- ipSymLS(R,w=W,eps=1e-15)

## ----ipsymlsplot, fig.cap = "WALS biplot of the correlation matrix of the wheat kernel data."----
bplot(Fp.als,Fp.als,rowch=NA,colch=NA,collab=colnames(R),
      xl=c(-1.1,1.1),yl=c(-1.1,1.1),main="WALS")
circle()

## -----------------------------------------------------------------------------
Rhat.wals <- Fp.als%*%t(Fp.als)
sqrt(diag(Rhat.wals))
rmse.als <- rmse(R,Rhat.wals)
rmse.als

## -----------------------------------------------------------------------------
p <- nrow(R)
W <- matrix(1,p,p)
diag(W) <- 0
out.wals <- wAddPCA(R, W, add = "nul", verboseout = FALSE, epsout = 1e-10)
Rhat.wals <- out.wals$a%*%t(out.wals$b)
out.eig <- eigen(Rhat.wals)
Fp.adj <- out.eig$vectors[,1:2]%*%diag(sqrt(out.eig$values[1:2]))
rmse.als <- rmse(R,Rhat.wals)
rmse.als

## -----------------------------------------------------------------------------
p <- nrow(R)
W <- matrix(1,p,p)
diag(W) <- 0
out.wals <- wAddPCA(R, W, add = "one", verboseout = FALSE, epsout = 1e-10)
delta <- out.wals$delta[1,1]
Rhat <- out.wals$a%*%t(out.wals$b)
out.eig <- eigen(Rhat)
Fp.adj <- out.eig$vectors[,1:2]%*%diag(sqrt(out.eig$values[1:2]))

## ----walsbiplot, fig.cap = "WALS biplot of the correlation matrix of the wheat kernel data, with the use of an additive adjustment."----
bplot(Fp.adj,Fp.adj,rowch=NA,colch=NA,collab=colnames(R),
      xl=c(-1.2,1.2),yl=c(-1.2,1.2),main="WALS adjusted")
circle()

## -----------------------------------------------------------------------------
Rhat.adj <- Fp.adj%*%t(Fp.adj)+delta
rmse.adj <- rmse(R,Rhat.adj)
rmse.adj

## -----------------------------------------------------------------------------
rmsevector <- c(rmse.crg,rmse.lin,rmse.pca,rmse.mds,rmse.pfa,rmse.als,rmse.adj)
methods <- c("CRG (cos)","CRG (lin)","PCA","MDS",
"PFA","WALS R","WALS Radj")
results <- data.frame(methods,rmsevector)
results <- results[order(rmsevector),]
results

## -----------------------------------------------------------------------------
output <- FitRwithPCAandWALS(R,eps=1e-15)
rmsePCAandWALS(R,output)

