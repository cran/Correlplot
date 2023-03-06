\name{FitRwithPCAandWALS}
\alias{FitRwithPCAandWALS}
\title{
  Calculate a low-rank approximation to the correlation matrix with four methods
}
\description{
  Function \code{FitRwithPCAandWALS} uses principal component analysis (PCA) and weighted alternating least squares (WALS) to
  calculate different low-rank approximations to the correlation matrix.
}
\usage{
FitRwithPCAandWALS(R, nd = 2, itmaxout = 10000, itmaxin = 10000, eps = 1e-08)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{R}{The correlation matrix}
  \item{nd}{The dimensionality of the low-rank solution (2 by default)}
  \item{itmaxout}{Maximum number of iterations for the outer loop of the algorithm}
  \item{itmaxin}{Maximum number of iterations for the inner loop of the algorithm}
  \item{eps}{Numerical criterion for convergence of the outer loop}
}
\details{
  Four methods are run succesively: standard PCA; PCA with an additive adjustment; WALS avoiding the fit of the diagonal;
  WALS avoiding the fit of the diagonal and with an additive adjustment.
}
\value{A list object with fields:
\item{Rhat.pca}{Low-rank approximation obtained by PCA}
\item{Rhat.pca.adj}{Low-rank approximation obtained by PCA with adjustment}
\item{Rhat.wals}{Low-rank approximation obtained by WALS without fitting the diagonal}
\item{Rhat.wals.adj}{Low-rank approximation obtained by WALS without fitting the diagonal and with adjustment}
}
\references{Graffelman, J. and De Leeuw, J. (2023) Improved approximation and visualization of the correlation matrix. The American Statistician pp. 1--20. Available online as latest article  \doi{10.1080/00031305.2023.2186952}}
\author{ Jan Graffelman (jan.graffelman@upc.edu) }
\seealso{
\code{\link{wAddPCA}}
}
\examples{
data(HeartAttack)
X <- HeartAttack[,1:7]
X[,7] <- log(X[,7])
colnames(X)[7] <- "logPR"
R <- cor(X)
\dontrun{
out <- FitRwithPCAandWALS(R)
}
}
\keyword{optimize}
\keyword{multivariate}
