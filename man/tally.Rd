\name{tally}
\alias{tally}
\title{
  Create a tally on a biplot vector
}
\description{
  Function \code{tally} marks of a set of dots on a biplot vector. It is thought for biplot vectors representing correlations,
  such that their correlation scale becomes visible, without doing a full calibration with tick marks and tick mark labels.
}
\usage{
tally(G, adj = 0, values = seq(-1, 1, by = 0.2), pch = 19, dotcolor = "black", cex = 0.5,
      color.negative = "red", color.positive = "blue")
}
\arguments{
  \item{G}{Matrix with biplot coordinates of the variables}
  \item{adj}{A scalar adjustment for the correlations}
  \item{values}{The values of the correlations to be marked off by dots}
  \item{pch}{The character code used for marking off correlations}
  \item{dotcolor}{The colour of the dots that are marked off}
  \item{cex}{The character expansion factor for a dot.}
  \item{color.negative}{The colour of the segments of the negative part of the correlation scale}
  \item{color.positive}{The colour of the segments of the positive part of the correlation scale}
}
\value{NULL}
\references{Graffelman, J. and De Leeuw, J. (2023) Improved approximation and visualization of the correlation matrix.
The American Statistician pp. 1--20. \doi{10.1080/00031305.2023.2186952}}
\author{ Jan Graffelman (jan.graffelman@upc.edu) }
\seealso{
    \code{\link{bplot}}, \code{\link{calibrate}} 
}
\examples{
data(goblets)
R <- cor(goblets)
results <- eigen(R)
V  <- results$vectors
Dl <- diag(results$values)
#
# Calculate correlation biplot coordinates
#
G  <- crossprod(t(V[,1:2]),sqrt(Dl[1:2,1:2]))
#
# Make the biplot
#
bplot(G,G,rowch=NA,colch=NA,collab=colnames(R),
      xl=c(-1.1,1.1),yl=c(-1.1,1.1))
#
# Create a correlation tally stick for variable X1
#
tally(G[1,])
}
\keyword{hplot}
\keyword{aplot}

