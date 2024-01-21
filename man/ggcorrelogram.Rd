\name{ggcorrelogram}
\alias{ggcorrelogram}
\title{
  Create a correlogram as a ggplot object.
}
\description{
  Function \code{ggcorrelogram} creates a correlogram of a correlation matrix using ggplot graphics.
}
\usage{
ggcorrelogram(R, labs = colnames(R), ifun = "cos", cex = 1, main = "", ntrials = 50,
              xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2), hjust = 1, vjust = 2, size = 2,
	      main.size = 8)
}
\arguments{
  \item{R}{a correlation matrix}
  \item{labs}{a vector of labels for the variables}
  \item{ifun}{the interpretation function ("cos" or "lincos")}
  \item{cex}{character expansion factor for the variable labels}
  \item{main}{a title for the correlogram}
  \item{ntrials}{number of starting points for the optimization routine}
  \item{xlim}{limits for the x axis (e.g. c(-1.2,1.2))}
  \item{ylim}{limits for the y axis (e.g. c(-1.2,1.2))}
  \item{hjust}{horizontal adjustment of variable labels (by default 1 for all variables)}
  \item{vjust}{vertical adjustment of variable labels (by default 2 for all variables)}
  \item{size}{font size for the labels of the variables}
  \item{main.size}{font size of the main title of the correlogram}
}
\details{
\code{ggcorrelogram} makes a correlogram on the basis of a set of
angles. All angles are given w.r.t the positive x-axis. Variables are
represented by unit vectors emanating from the origin.
}
\value{
  A ggplot object. Field \code{theta} of the output contains the angles for the variables. 
}
\references{Trosset, M.W. (2005) Visualizing correlation. Journal of Computational and Graphical Statistics
  14(1), pp. 1--19}
\author{ Jan Graffelman (jan.graffelman@upc.edu) }
\seealso{\code{\link{correlogram}},\code{\link{fit_angles}},\code{\link{nlminb}}}
\examples{
 set.seed(123)
 X <- matrix(rnorm(90),ncol=3)
 R <- cor(X)
 angles <- ggcorrelogram(R)
}
\keyword{hplot}

