\name{ggtally}
\alias{ggtally}
\title{
  Create a correlation tally stick on a biplot vector
}
\description{
  Function \code{ggtally} puts a series of dots along a biplot vector of a correlation matrix,
  so marking the change in correlation along the vector with specified values.
}
\usage{
ggtally(G, p1, adj = 0, values = seq(-1, 1, by = 0.2), dotsize = 0.1, dotcolour = "black")
}
\arguments{
  \item{G}{A matrix (or vector) of biplot markers}
  \item{p1}{A ggplot2 object with a biplot}
  \item{adj}{A scalar adjustment for the correlations}
  \item{values}{Values of the correlations to be marked off by dots}
  \item{dotsize}{Size of the dot}
  \item{dotcolour}{Colour of the dot}
}
\details{
  Any set of values for the correlation to be marked off can be used, though a standard scale with 0.2 increments is recommmended. 
}
\value{A ggplot2 object with the updated biplot}
\references{
Graffelman, J. and De Leeuw, J. (2023) On the visualisation of the correlation matrix. Available online. \doi{10.48550/arXiv.2211.13150}}
\author{ Jan Graffelman (jan.graffelman@upc.edu) }
\seealso{
   \code{\link{ggbplot}}
}
\examples{
library(calibrate)
data(goblets)
R <- cor(goblets)
out.sd <- eigen(R)
V  <- out.sd$vectors[,1:2]
Dl <- diag(out.sd$values[1:2])
Gp <- crossprod(t(V),sqrt(Dl))
pca.df <- data.frame(Gp)
pca.df$strings <- colnames(R)
colnames(pca.df) <- c("PA1","PA2","strings")
p1 <- ggbplot(pca.df,pca.df,main="PCA correlation biplot",xlab="",ylab="",rowarrow=TRUE,
              rowcolor="blue",rowch="",colch="")
p1 <- ggtally(Gp,p1,values=seq(-0.2,0.6,by=0.2),dotsize=0.1)
}
\keyword{hplot}

