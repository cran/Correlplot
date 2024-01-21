`jointlim` <-
function(X,Y) {
   xma <- max(c(X[,1],Y[,1]))
   xmi <- min(c(X[,1],Y[,1]))
   yma <- max(c(X[,2],Y[,2]))
   ymi <- min(c(X[,2],Y[,2]))
   xlim <- c(xmi,xma)
   ylim <- c(ymi,yma)
   return(list(xlim=xlim,ylim=ylim))
}

