\name{lincos}
\alias{lincos}
\title{
Linearized cosine function
}
\description{
Function \code{lincos} linearizes the cosine function over the interval
[0,2pi]. The function returns -2/pi*x + 1 over [0,pi] and 2/pi*x - 3
over [pi,2pi]  
}
\usage{
lincos(x)
}
\arguments{
  \item{x}{angle in radians}
}
\value{a real number in [-1,1].}
\references{
  Graffelman, J. (2012) Linear-angle correlation plots: new graphs for revealing correlation structure. Journal of Computational and Graphical Statistics. 22(1): 92-106. 
}
\author{Jan Graffelman (jan.graffelman@upc.edu)}
\seealso{\code{\link{cos}}}
\examples{
angle <- pi
y <- lincos(angle)
print(y)
}
\keyword{misc}

