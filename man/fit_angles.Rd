\name{fit_angles}
\alias{fit_angles}
\title{Fit angles to a correlation matrix}
\description{
\code{fit_angles} finds a set of optimal angles for representing a
particular correlation matrix by angles between vectors
}
\usage{
fit_angles(R, ifun = "cos", ntrials = 10, verbose = FALSE)
}
\arguments{
  \item{R}{a correlation matrix.}
  \item{ifun}{an angle interpretation function (cosine, by default).}
  \item{ntrials}{number of trials for optimization routine \code{nlminb}}
  \item{verbose}{be silent (FALSE), or produce more output (TRUE)}
}
\value{a vector of angles (in radians)}
\references{Trosset, M.W. (2005) Visualizing correlation. Journal of Computational and Graphical Statistics
14(1), pp. 1--19}
\author{ anonymous }
\seealso{\code{\link{nlminb}}}
\examples{
X <- matrix(rnorm(90),ncol=3)
R <- cor(X)
angles <- fit_angles(R)
print(angles)
}
\keyword{optimize}
