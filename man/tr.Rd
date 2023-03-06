\name{tr}
\alias{tr}
\title{Compute the trace of a matrix}
\description{
   \code{tr} computes the trace of a matrix.
}
\usage{
tr(X)
}
\arguments{
  \item{X}{a (square) matrix}
}
\value{
  the trace (a scalar)
}
\author{ Jan Graffelman (jan.graffelman@upc.edu) }
\examples{
X <- matrix(runif(25),ncol=5)
print(X)
print(tr(X))
}
\keyword{array}
\keyword{algebra}

