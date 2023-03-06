\name{rmsePCAandWALS}
\alias{rmsePCAandWALS}
\title{
  Generate a table of root mean square error (RMSE) statistics for principal component analysis (PCA) and
  weighted alternating least squares (WALS).
}
\description{
  Function \code{rmsePCAandWALS} creates table with the RMSE for each variable, for a low-rank
  approximation to the correlation matrix obtained by PCA or WALS.
}
\usage{
rmsePCAandWALS(R, output, digits = 4, diagonals = c(FALSE,FALSE,TRUE,TRUE))
}
\arguments{
  \item{R}{The correlation matrix}
  \item{output}{A list object with four approximationst to the correlation matrix}
  \item{digits}{The number of digits used in the output}
  \item{diagonals}{Vector of four logicals for omitting the diagonal of the correlation matrix
  for RMSE calculations. Defaults to c(FALSE,FALSE,TRUE,TRUE), to include the diagonal
  for PCA and exclude it for WALS}
}
\value{A matrix with one row per variable and four columns for RMSE statistics.}
\references{
  Graffelman, J. and De Leeuw, J. (2023) Improved approximation and visualization of the correlation matrix. The American Statistician pp. 1--20. Available online as latest article  \doi{10.1080/00031305.2023.2186952}
}
\author{ Jan Graffelman (jan.graffelman@upc.edu) }
\seealso{
  \code{\link{FitRwithPCAandWALS}}
}
\examples{
data(HeartAttack)
X <- HeartAttack[,1:7]
X[,7] <- log(X[,7])
colnames(X)[7] <- "logPR"
R <- cor(X)
\dontrun{
out <- FitRwithPCAandWALS(R)
Results <- rmsePCAandWALS(R,out)
}
}
\keyword{optimize}
\keyword{multivariate}

