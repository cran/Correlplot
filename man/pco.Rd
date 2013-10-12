\name{pco}
\alias{pco}
\title{Principal Coordinate Analysis}
\description{
  \code{pco} is a program for Principal Coordinate Analysis.
}
\usage{
pco(Dis)
}
\arguments{
  \item{Dis}{A distance or dissimilarity matrix}
}
\details{
The program \code{pco} does a principal coordinates analysis of a 
dissimilarity (or distance) matrix (Dij) where the diagonal elements, 
Dii, are zero.

Note that when we dispose of a similarity matrix rather that a distance
matrix, a transformation is needed before calling coorprincipal. For 
instance, if Sij is a similarity matrix, Dij might be obtained as
Dij = 1 - Sij/diag(Sij)

Goodness of fit calculations need to be revised such as to deal
(in different ways) with negative eigenvalues.
}
\value{
  \item{PC}{the principal coordinates}
  \item{Dl}{all eigenvalues of the solution}
  \item{Dk}{the positive eigenvalues of the solution}
  \item{B}{double centred matrix for the eigenvalue decomposition}
  \item{decom}{the goodness of fit table}
}
\author{ Jan Graffelman (jan.graffelman@upc.edu) }
\seealso{ \code{\link{cmdscale}} }
\examples{
citynames <- c("Aberystwyth","Brighton","Carlisle","Dover","Exeter","Glasgow","Hull",
"Inverness","Leeds","London","Newcastle", "Norwich")    
A <-matrix(c(
0,244,218,284,197,312,215,469,166,212,253,270,
244,0,350,77,167,444,221,583,242,53,325,168,
218,350,0,369,347,94,150,251,116,298,57,284,
284,77,369,0,242,463,236,598,257,72,340,164,
197,167,347,242,0,441,279,598,269,170,359,277,
312,444,94,463,441,0,245,169,210,392,143,378,
215,221,150,236,279,245,0,380,55,168,117,143,
469,583,251,598,598,169,380,0,349,531,264,514,
166,242,116,257,269,210,55,349,0,190,91,173,
212,53,298,72,170,392,168,531,190,0,273,111,
253,325,57,340,359,143,117,264,91,273,0,256,
270,168,284,164,277,378,143,514,173,111,256,0),ncol=12)
rownames(A) <- citynames
colnames(A) <- citynames
out <- pco(A)
plot(out$PC[,2],-out$PC[,1],pch=19,asp=1)
textxy(out$PC[,2],-out$PC[,1],rownames(A))
}
\keyword{multivariate}

