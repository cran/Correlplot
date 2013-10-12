\name{linangplot}
\alias{linangplot}
\title{Linang plot}
\description{
   \code{linangplot} produces a plot of two variables, such that the correlation between the two variables is linear in the angle.
}
\usage{
linangplot(x, y, tmx = NULL, tmy = NULL, ...)
}
\arguments{
  \item{x}{x variable}
  \item{y}{y variable}
  \item{tmx}{vector of tickmarks for the x variable}
  \item{tmy}{vector of tickmarks for the y variable}
  \item{\dots}{additional arguments for the plot routine}
}
\value{
  \item{Xt}{coordinates of the points}
  \item{B}{axes for the plot}
  \item{r}{correlation coefficient}
  \item{angledegrees}{angle between axes in degrees}
  \item{angleradians}{angle between axes in radians}
  \item{r}{correlation coefficient}
}
\author{ Jan Graffelman (jan.graffelman@upc.edu) }
\seealso{
\code{\link{plot}}\code{\link{correlogram}}
}
\examples{
x <- runif(10)
y <- rnorm(10)
linangplot(x,y)
}
\keyword{hplot}

