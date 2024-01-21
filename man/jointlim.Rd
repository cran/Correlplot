\name{jointlim}
\alias{jointlim}
\title{Establish limits for x and y axis}
\description{
  jointlim computes a sensible range for x and y axis if two sets of
  points are to be plotted simultaneously
}
\usage{
jointlim(X, Y)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{Matrix of coordinates}
  \item{Y}{Matrix of coordinates}
}
\value{
  \item{xlim}{minimum and maximum for x-range}
  \item{ylim}{minimum and maximum for y-range}
}
\author{ Jan Graffelman (jan.graffelman@upc.edu) }
\examples{
X <- matrix(runif(20),ncol=2)
Y <- matrix(runif(20),ncol=2)
print(jointlim(X,Y)$xlim)
}
\keyword{misc}


