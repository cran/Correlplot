## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
knitr::opts_chunk$set(fig.align = "center")
knitr::opts_chunk$set(fig.width = 6, fig.height = 6) 

## ----preinstall---------------------------------------------------------------
library(calibrate)
library(ggplot2)
library(corrplot)
library(Correlplot)

## -----------------------------------------------------------------------------
data("Kernels")
X <- Kernels[Kernels$variety==1,]
X <- X[,-8]
head(X)

## -----------------------------------------------------------------------------
p <- ncol(X)
R <- cor(X)
round(R,digits=3)

## ----corrgram, fig.cap = "A corrgram of the wheat kernel data."---------------
corrplot(R, method="circle",type="lower")

## ----pnames, echo = FALSE-----------------------------------------------------
vnames <- substr(colnames(R),1,4)
colnames(R) <- rownames(R) <- vnames
colnames(X) <- vnames

## ----correlogram, fig.cap = "The correlogram of the wheat kernel data."-------
theta.cos <- correlogram(R,xlim=c(-1.1,1.1),ylim=c(-1.1,1.1),main="CRG")

## -----------------------------------------------------------------------------
Rhat.cor <- angleToR(theta.cos)

## -----------------------------------------------------------------------------
W1 <- matrix(1,p,p)
rmse.crg <- rmse(R,Rhat.cor,W=W1)
rmse.crg

## ----linearcorrelogram, fig.cap="Linear correlogram of the wheat kernel data."----
theta.lin <- correlogram(R,ifun="lincos",labs=colnames(R),xlim=c(-1.1,1.1),
                         ylim=c(-1.1,1.1),main="CRG")

## -----------------------------------------------------------------------------
Rhat.corlin <- angleToR(theta.lin,ifun="lincos")
rmse.lin <- rmse(R,Rhat.corlin,W=W1)
rmse.lin

## ----ggcorrelogram, fig.cap="Correlogram of the wheat kernel data on a ggplot2 canvas."----
set.seed(123)
crgR <- ggcorrelogram(R,main="CRG",vjust=c(0,2,-1,2,0,-1,2),
                     hjust=0)
crgR
crgR$theta

## ----pcabiplot, fig.cap = "PCA biplot of the (standardized) wheat kernel data."----
n <- nrow(X)
Xt <- scale(X)/sqrt(n-1)
res.svd <- svd(Xt)
Fs <- sqrt(n)*res.svd$u # standardized principal components
Gp <- res.svd$v%*%diag(res.svd$d) # biplot coordinates for variables
bplot(Fs,Gp,colch=NA,collab=colnames(X),xlab="PC1",ylab="PC2",main="PCA")

## ----talliedcorrelationcircle, fig.cap = "PCA biplot of the correlation matrix with correlation tally sticks."----
bplot(Gp,Gp,colch=NA,rowch=NA,collab=colnames(X),xl=c(-1,1),yl=c(-1,1),main="PCA")
circle()
tally(Gp[-6,1:2],values=seq(-0.2,0.8,by=0.2))

## -----------------------------------------------------------------------------
Rhat.pca <- Gp[,1:2]%*%t(Gp[,1:2])

## -----------------------------------------------------------------------------
rmse.pca <- rmse(R,Rhat.pca,W=W1)
rmse.pca

## -----------------------------------------------------------------------------
rmse(R,Rhat.pca,W=W1,per.variable=TRUE)

## -----------------------------------------------------------------------------
limits <- jointlim(Fs,Gp)
limits

## -----------------------------------------------------------------------------
df.rows <- data.frame(Fs[,1:2])
colnames(df.rows) <- c("PA1","PA2")
df.rows$strings <- 1:n

df.cols <- data.frame(Gp[,1:2])
colnames(df.cols) <- c("PA1","PA2")
df.cols$strings <- colnames(R)

## -----------------------------------------------------------------------------
lambda      <- res.svd$d^2
lambda.frac <- lambda/sum(lambda)
lambda.cum  <- cumsum(lambda.frac)
decom <- round(rbind(lambda,lambda.frac,lambda.cum),digits=3)
colnames(decom) <- paste("PC",1:p,sep="")
decom

## -----------------------------------------------------------------------------
xlab <- paste("PC1 (",round(100*lambda.frac[1],digits=1),"%)",sep="")
ylab <- paste("PC2 (",round(100*lambda.frac[2],digits=1),"%)",sep="")

## ----ggbiplot, fig.cap = "PCA biplot on a ggplot canvas"----------------------
biplotX <- ggbplot(df.rows,df.cols,main="PCA biplot",xlab=xlab,
              ylab=ylab,xlim=limits$xlim,ylim=limits$ylim,
              colch="")
biplotX

## -----------------------------------------------------------------------------
lambdasq      <- lambda^2
lambdasq.frac <- lambdasq/sum(lambdasq)
lambdasq.cum  <- cumsum(lambdasq.frac)
decomsq <- round(rbind(lambdasq,lambdasq.frac,lambdasq.cum),
                 digits=3)
colnames(decomsq) <- paste("PC",1:p,sep="")
decomsq

xlab <- paste("PC1 (",round(100*lambdasq.frac[1],digits=1),"%)",sep="")
ylab <- paste("PC2 (",round(100*lambdasq.frac[2],digits=1),"%)",sep="")


## ----anotherbiplot, fig.cap = "PCA Correlation biplot on a ggplot canvas."----
biplotR <- ggbplot(df.cols,df.cols,main="PCA Correlation biplot",xlab=xlab,
              ylab=ylab,xlim=c(-1,1),ylim=c(-1,1),
              rowarrow=TRUE,rowcolor="blue",colch="",
              rowch="")

biplotR

## ----yetanotherbiplot, fig.cap = "PCA Correlation biplot with correlation tally sticks."----
biplotR <- ggtally(Gp[-6,1:2],biplotR,values=seq(-0.2,0.8,by=0.2),dotsize=0.2)
biplotR <- ggtally(Gp[6,1:2],biplotR,values=seq(-0.01,0.01,by=0.01),dotsize=0.2)
biplotR

## ----mdsplot, fig.cap = "MDS map of the correlation matrix of the wheat kernel data."----
Di <- sqrt(2*(1-R))
out.mds <- cmdscale(Di,eig = TRUE)
Fp <- out.mds$points[,1:2]
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
rmse.mds <- rmse(R,Rhat.mds,W=W1)
rmse.mds

## ----mdsggplot, fig.cap = "MDS map on ggplot canvas."-------------------------
colnames(Fp) <- paste("PA",1:2,sep="")
rownames(Fp) <- as.character(1:nrow(Fp))
Fp <- data.frame(Fp)
Fp$strings <- colnames(R)
MDSmap <- ggplot(Fp, aes(x = PA1, y = PA2)) + 
            coord_equal(xlim = c(-1,1), ylim = c(-1,1)) +
            ggtitle("MDS") +
            xlab("First principal axis") + ylab("Second principal axis") +
            theme(plot.title = element_text(hjust = 0.5, 
                                            size = 8),
                  axis.ticks = element_blank(),
                  axis.text.x = element_blank(),
                  axis.text.y = element_blank())
MDSmap <- MDSmap + geom_point(data = Fp, aes(x = PA1, y = PA2), 
                        colour = "black", shape = 1)
MDSmap <- MDSmap + geom_text(data = Fp, aes(label = strings), 
                        size = 3, alpha = 1,
                      vjust = 2) 

Z <- matrix(NA,nrow=nrow(ii),ncol=4)
for(i in 1:nrow(ii)) {
  Z[i,] <- c(Fp[ii[i,1],1],Fp[ii[i,1],2],Fp[ii[i,2],1],Fp[ii[i,2],2])
}
colnames(Z) <- c("X1","Y1","X2","Y2")
Z <- data.frame(Z)

MDSmap <- MDSmap + geom_segment(data=Z,aes(x=X1,y=Y1,
                                           xend=X2,yend=Y2),
                         alpha=0.75,color="red",linetype=2)   
MDSmap

## ----pfa----------------------------------------------------------------------
out.pfa <- Correlplot::pfa(X,verbose=FALSE)
L <- out.pfa$La

## -----------------------------------------------------------------------------
Rhat.pfa <- L[,1:2]%*%t(L[,1:2])
rmse.pfa <- rmse(R,Rhat.pfa)
rmse.pfa

## ----pfabiplot, fig.cap = "PFA biplot of the correlation matrix of the wheat kernel data."----
bplot(L,L,pch=NA,xl=c(-1,1),yl=c(-1,1),xlab="Factor 1",ylab="Factor 2",main="PFA",rowch=NA,
      colch=NA)
circle()

## -----------------------------------------------------------------------------
diag(out.pfa$Psi)

## ----pfaggbiplot, fig.cap = "PFA biplot on a ggplot canvas."------------------
L.df <- data.frame(L,rownames(L))
colnames(L.df) <- c("PA1","PA2","strings")
ggbplot(L.df,L.df,xlab="Factor 1",ylab="Factor 2",main="PFA biplot",
        rowarrow=TRUE,rowcolor="blue",colch="",rowch="")

## -----------------------------------------------------------------------------
Wdiag0 <- matrix(1,nrow(R),nrow(R))
diag(Wdiag0) <- 0
Fp.als <- ipSymLS(R,w=Wdiag0,eps=1e-15)

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
out.wals <- wAddPCA(R, Wdiag0, add = "nul", verboseout = FALSE, epsout = 1e-10)
Rhat.wals <- out.wals$a%*%t(out.wals$b)
out.eig <- eigen(Rhat.wals)
Fp.adj <- out.eig$vectors[,1:2]%*%diag(sqrt(out.eig$values[1:2]))
rmse.als <- rmse(R,Rhat.wals)
rmse.als

## ----ipsymlsggplot, fig.cap = "WALS biplot of the correlation matrix of the wheat kernel data."----
Fp.als.df <- data.frame(Fp.als,colnames(R))
colnames(Fp.als.df) <- c("PA1","PA2","strings")
ggbplot(Fp.als.df,Fp.als.df,xlab="Dimension 1",ylab="Dimension 2",main="WALS biplot",
        rowarrow=TRUE,rowcolor="blue",colch="",rowch="")

## -----------------------------------------------------------------------------
out.wals <- wAddPCA(R, Wdiag0, add = "one", verboseout = FALSE, epsout = 1e-10)
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

## ----walsggbiplot, fig.cap = "WALS biplot of the correlation matrix of the wheat kernel data, with additive adjustment and tally sticks."----
Fp.adj.df <- data.frame(Fp.adj,colnames(R))
colnames(Fp.adj) <- c("PA1","PA2")
colnames(Fp.adj.df) <- c("PA1","PA2","strings")
WALSbiplot.adj <- ggbplot(Fp.adj.df,Fp.adj.df,xlab="Dimension 1",ylab="Dimension 2",
                          main="WALS adjusted",rowarrow=TRUE,
                          rowcolor = "blue",rowch="",colch="")
WALSbiplot.adj <- ggtally(Fp.adj[-6,1:2],WALSbiplot.adj,
                          adj=-out.wals$delta[1,1],
                   values=seq(-0.2,0.8,by=0.2),dotsize=0.2)
WALSbiplot.adj

## -----------------------------------------------------------------------------
out.walsq.sym <- FitRDeltaQSym(R,Wdiag0,itmax.inner=30000,itmax.outer=30000,
                     eps=1e-8)
Gq.sym <- out.walsq.sym$C
rownames(Gq.sym) <- colnames(R)
Rhat.wsym <- out.walsq.sym$Rhat
rmse.wsym <- rmse(R,Rhat.wsym)
rmse.wsym

## ----newplot, fig.cap = "WALS biplot of the correlation matrix with column adjustment and tally sticks."----
Gq.sym.df <- data.frame(Gq.sym)
Gq.sym.df$strings <- colnames(R)
colnames(Gq.sym.df) <- c("PA1","PA2","strings")
ll <- 1.5
WALSbiplotq.sym <- ggbplot(Gq.sym.df,Gq.sym.df,main="WALS-q-sym biplot",xlab="Dimension 1",
                           ylab="Dimension 2",
                           ylim=c(-ll,ll),xlim=c(-ll,ll),circle=TRUE,
                           rowarrow=TRUE,rowcolor="blue",rowch="",
                           colch="")
for(i in c(1:5,7)) {
   WALSbiplotq.sym <- ggtally(Gq.sym[i,1:2],WALSbiplotq.sym,
                              adj=-out.walsq.sym$delta-out.walsq.sym$q[i], 
              values=seq(-0.2,1,by=0.2),dotsize=0.2)  
}
WALSbiplotq.sym

## -----------------------------------------------------------------------------
out.walsq.sym$delta + out.walsq.sym$q

