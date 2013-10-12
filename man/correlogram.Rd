\name{correlogram}
\alias{correlogram}
\title{
Plot a correlogram
}
\description{
\code{correlogram} plots a correlogram for a correlation matrix.
}
\usage{
correlogram(R,labs=colnames(R),ifun="cos",cex=1,main="",ntrials=10,
            xlim=c(-1.2,1.2),ylim=c(-1.2,1.2),...)
}
\arguments{
  \item{R}{a correlation matrix.}
  \item{labs}{a vector of labels for the variables.}
  \item{ifun}{the interpretation function ("cos" or "lincos")}
  \item{cex}{character expansion factor for the variable labels}
  \item{main}{a title for the correlogram}
  \item{ntrials}{number of starting points for the optimization routine}
  \item{xlim}{limits for the x axis (e.g. c(-1.2,1.2))}
  \item{ylim}{limits for the y axis (e.g. c(-1.2,1.2))}
  \item{\dots}{additional arguments for the \code{plot} function.}
}
\details{
\code{correlogram} makes a correlogram on the basis of a set of
angles. All angles are given w.r.t the positive x-axis. Variables are
represented by unit vectors emanating from the origin.
}
\value{NULL}
\references{Trosset, M.W. (2005) Visualizing correlation. Journal of Computational and Graphical Statistics
  14(1), pp. 1--19}
\author{ anonymous }
\seealso{
\code{\link{fit_angles}},\code{\link{nlminb}}
}
\examples{
X <- matrix(rnorm(90),ncol=3)
R <- cor(X)
correlogram(R)
}
\keyword{hplot}

