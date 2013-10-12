\name{angleToR}
\alias{angleToR}
\title{
Convert angles to correlations.
}
\description{
Function \code{angleToR} converts a vector of angles (in radians) to an
estimate of the correlation matrix, given an interpretation function.
}
\usage{
angleToR(x, ifun = "cos")
}
\arguments{
  \item{x}{a vector of angles (in radians)}
  \item{ifun}{the interpretation function ("cos" or "lincos")}
}
\value{A correlation matrix}
\references{
  Graffelman, J. (2012) Linear-angle correlation plots: new graphs for revealing correlation structure. Journal of Computational and Graphical Statistics. 22(1): 92-106. 
}
\author{Jan Graffelman (jan.graffelman@upc.edu)}
\seealso{\code{\link{cos}},\code{\link{lincos}}}
\examples{
angles <- c(0,pi/3)
R <- angleToR(angles)
print(R)
}
\keyword{misc}
