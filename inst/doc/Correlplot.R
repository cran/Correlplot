### R code from vignette source 'Correlplot.Rnw'

###################################################
### code chunk number 1: Correlplot.Rnw:78-80
###################################################
#install.packages("Correlplot")
library("Correlplot")


###################################################
### code chunk number 2: Correlplot.Rnw:86-92
###################################################
#install.packages("calibrate")
#install.packages("corrplot")
#install.packages("xtable")
library(calibrate)
library(corrplot)
library(xtable)


###################################################
### code chunk number 3: Correlplot.Rnw:107-112
###################################################
library(Correlplot)
data("Kernels")
X <- Kernels[Kernels$variety==1,]
X <- X[,-8]
head(X)


###################################################
### code chunk number 4: Correlplot.Rnw:117-119
###################################################
R <- cor(X)
xtable(R,digits=3)


###################################################
### code chunk number 5: Correlplot.Rnw:128-132
###################################################
#install.packages("corrplot")
library(corrplot)
R <- cor(X)
corrplot(R, method="circle",type="lower")


###################################################
### code chunk number 6: Correlplot.Rnw:144-146
###################################################
theta.cos <- correlogram(R,main="Correlogram wheat kernels",
            xlim=c(-1.3,1.3),ylim=c(-1.3,1.3))


###################################################
### code chunk number 7: Correlplot.Rnw:151-152
###################################################
Rhat.cor <- angleToR(theta.cos)


###################################################
### code chunk number 8: Correlplot.Rnw:157-158
###################################################
rmse.crg <- rmse(R,Rhat.cor,verbose=TRUE)


###################################################
### code chunk number 9: Correlplot.Rnw:167-170
###################################################
theta.lin <- correlogram(R,ifun="lincos",labs=colnames(R),
                         main="Linear Correlogram",
                         xlim=c(-1.3,1.3),ylim=c(-1.3,1.3))


###################################################
### code chunk number 10: Correlplot.Rnw:177-179
###################################################
Rhat.corlin <- angleToR(theta.lin,ifun="lincos")
rmse.lin <- rmse(R,Rhat.corlin,verbose=TRUE)


###################################################
### code chunk number 11: biplot
###################################################
n <- nrow(X)
Xt <- scale(X)/sqrt(n)
res.svd <- svd(Xt)
Fs <- sqrt(n)*res.svd$u
Gp <- res.svd$v%*%diag(res.svd$d)
bplot(Fs,Gp,colch=NA,collab=colnames(X),
      xlab = "First principal component",
      ylab="Second principal component")


###################################################
### code chunk number 12: Correlplot.Rnw:210-213
###################################################
bplot(Gp,Gp,colch=NA,rowch=NA,collab=colnames(X),
      xl=c(-1,1),yl=c(-1,1))
circle()


###################################################
### code chunk number 13: Correlplot.Rnw:219-221
###################################################
Rhat.pca <- Gp[,1:2]%*%t(Gp[,1:2])
rmse.pca <- rmse(R,Rhat.pca,verbose=TRUE)


###################################################
### code chunk number 14: Correlplot.Rnw:235-248
###################################################
Di <- sqrt(2*(1-R))
out.mds <- cmdscale(Di,eig = TRUE)
Fp <- out.mds$points
plot(Fp[,1],Fp[,2],asp=1,xlab="First principal axis",
     ylab="Second principal axis")
textxy(Fp[,1],Fp[,2],colnames(R),cex=0.75)

ii <- which(R < 0,arr.ind = TRUE)

for(i in 1:nrow(ii)) {
  segments(Fp[ii[i,1],1],Fp[ii[i,1],2],
           Fp[ii[i,2],1],Fp[ii[i,2],2],col="red",lty="dashed")
}


###################################################
### code chunk number 15: Correlplot.Rnw:254-257
###################################################
Dest <- as.matrix(dist(Fp[,1:2]))
Rhat.mds <- 1-0.5*Dest*Dest
rmse.mds <- rmse(R,Rhat.mds,verbose=TRUE)


###################################################
### code chunk number 16: Correlplot.Rnw:266-268
###################################################
out.pfa <- pfa(X)
L <- out.pfa$La


###################################################
### code chunk number 17: Correlplot.Rnw:274-276
###################################################
Rhat.pfa <- L[,1:2]%*%t(L[,1:2])
rmse.pfa <- rmse(R,Rhat.pfa,verbose=TRUE)


###################################################
### code chunk number 18: pfaplot
###################################################
opar <- par(bty="n",xaxt="n",yaxt="n")
plot(L[,1],L[,2],pch=NA,asp=1,xlim=c(-1,1),ylim=c(-1,1),
     xl="Factor 1",yl="Factor 2")
origin()
arrows(0,0,L[,1],L[,2],angle=10,length=0.1,col="blue")
textxy(L[,1],L[,2],colnames(X),cex=1)
circle()
par(opar)


###################################################
### code chunk number 19: Correlplot.Rnw:304-307
###################################################
W <- matrix(1,nrow(R),nrow(R))
diag(W) <- 0
Fp.als <- ipSymLS(R,w=W,eps=1e-15)


###################################################
### code chunk number 20: Correlplot.Rnw:311-314
###################################################
bplot(Fp.als,Fp.als,rowch=NA,colch=NA,collab=colnames(R),
      xl=c(-1.2,1.2),yl=c(-1.2,1.2),main="WALS")
circle()


###################################################
### code chunk number 21: Correlplot.Rnw:319-321
###################################################
Rhat.wals <- Fp.als%*%t(Fp.als)
rmse.als <- rmse(R,Rhat.wals,verbose=TRUE)


###################################################
### code chunk number 22: Correlplot.Rnw:338-342
###################################################
delta <- 0.07
W <- matrix(1,nrow(R),nrow(R))
diag(W) <- 0
Fp.adj <- ipSymLS(R-delta,w=W,verbose=FALSE,eps=1e-10,itmax=1000)


###################################################
### code chunk number 23: Correlplot.Rnw:348-351
###################################################
bplot(Fp.adj,Fp.adj,rowch=NA,colch=NA,collab=colnames(R),
      xl=c(-1.3,1.3),yl=c(-1.3,1.3),main="WALS adjusted")
circle()


###################################################
### code chunk number 24: Correlplot.Rnw:357-359
###################################################
Rhat.adj <- Fp.adj%*%t(Fp.adj) + delta
rmse.adj <- rmse(R,Rhat.adj,verbose=TRUE)


###################################################
### code chunk number 25: Correlplot.Rnw:368-372
###################################################
rmsevector <- c(rmse.crg,rmse.lin,rmse.pca,rmse.mds,rmse.pfa,rmse.als,rmse.adj)
methods <- c("Correlogram (cosine)","Correlogram (linear)","PCA","MDS",
"PFA","WALS R","WALS Radj")
xtable(data.frame(methods,rmsevector),digits=c(0,0,4))


