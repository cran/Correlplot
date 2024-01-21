\name{ggbplot}
\alias{ggbplot}
\title{
  Create a biplot with ggplot2
}
\description{
  Function \code{ggbiplot} creates a biplot of a matrix with ggplot2 graphics.
}
\usage{
ggbplot(A, B, main = "", circle = TRUE, xlab = "", ylab = "", main.size = 8,
xlim = c(-1, 1), ylim = c(-1, 1), rowcolor = "red", rowch = 1, colcolor = "blue",
colch = 1, rowarrow = FALSE, colarrow = TRUE)
}
\arguments{
  \item{A}{A dataframe with coordinates and names for the biplot row markers}
  \item{B}{A dataframe with coordinates and names for the biplot column markers}
  \item{main}{A title for the biplot}
  \item{circle}{Draw a unit circle (\code{circle=TRUE}) or not (\code{circle=FALSE})}
  \item{xlab}{The label for the x axis}
  \item{ylab}{The label for the y axis}
  \item{main.size}{Size of the main title}
  \item{xlim}{Limits for the horizontal axis}
  \item{ylim}{Limits for the vertical axis}
  \item{rowcolor}{Color used for the row markers}
  \item{rowch}{Symbol used for the row markers}
  \item{colcolor}{Color used for the column markers}
  \item{colch}{Symbol used for the column markers}
  \item{rowarrow}{Draw arrows from the origin to the row markers (\code{rowarrow=TRUE}) or not}
  \item{colarrow}{Draw arrows from the origin to the column markers (\code{colarrow=TRUE}) or not}
}
\details{
  Dataframes \code{A} and \code{B} must consists of three columns labeled "PA1", "PA2" (coordinates of the first and second principal axis) and a column "strings" with the labels for the coordinates.
  
  Dataframe \code{B} is optional. If it is not specified, a biplot with a single set of markers is constructed, for which the row settings must be specified.}
\value{A ggplot2 object}
\references{
Graffelman, J. and De Leeuw, J. (2023) On the visualisation of the correlation matrix. Available online. \doi{10.48550/arXiv.2211.13150}}
\author{ Jan Graffelman (jan.graffelman@upc.edu) }
\seealso{
  \code{\link{bplot}},\code{\link{ggtally}},\code{\link{biplot}}
}
\examples{
data("HeartAttack")
X <- as.matrix(HeartAttack[,1:7])
n <- nrow(X)
Xt <- scale(X)/sqrt(n-1)
res.svd <- svd(Xt)
Fs <- sqrt(n)*res.svd$u # standardized principal components
Gp <- crossprod(t(res.svd$v),diag(res.svd$d)) # biplot coordinates for variables
rows.df <- data.frame(Fs[,1:2],as.character(1:n))
colnames(rows.df) <- c("PA1","PA2","strings")
cols.df <- data.frame(Gp[,1:2],colnames(X))
colnames(cols.df) <- c("PA1","PA2","strings")
ggbplot(rows.df,cols.df,xlab="PA1",ylab="PA2",main="PCA")
}
\keyword{hplot}

