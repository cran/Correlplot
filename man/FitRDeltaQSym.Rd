\name{FitRDeltaQSym}
\alias{FitRDeltaQSym}
\title{Approximation of a correlation matrix with column adjustment and symmetric low rank factorization}
\description{
Program \code{FitRDeltaQSym} calculates a low rank factorization for a correlation matrix. It adjusts for column effects, and the approximation is therefore asymmetric.  
}
\usage{
FitRDeltaQSym(R, W = NULL, nd = 2, eps = 1e-10, delta = 0, q = colMeans(R),
              itmax.inner = 1000, itmax.outer = 1000, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{R}{A correlation matrix}
  \item{W}{A weight matrix (optional)}
  \item{nd}{The rank of the low rank approximation}
  \item{eps}{The convergence criterion}
  \item{delta}{Initial value for the scalar adjustment (zero by default)}
  \item{q}{Initial values for the column adjustments (random by default)}
  \item{itmax.inner}{Maximum number of iterations for the inner loop of the algorithm}
  \item{itmax.outer}{Maximum number of iterations for the outer loop of the algorithm}
  \item{verbose}{Print information or not}
}
\details{
Program \code{FitRDeltaQSym} implements an iterative algorithm for the low rank factorization of the correlation matrix. It decomposes the correlation matrix as R = delta J + 1 q' + G G' + E. The approximation of R is ultimately asymmetric, but the low rank factorization used for biplotting (G G') is symmetric. 
}
\value{A list object with fields:
\item{delta}{The final scalar adjustment}
\item{Rhat}{The final approximation to the correlation matrix}
\item{C}{The matrix of biplot vectors}
\item{rmse}{The root mean squared error}
\item{q}{The final column adjustments}
}
\references{
Graffelman, J. and De Leeuw, J. (2023) Improved approximation and visualization of the correlation matrix. The American Statistician pp. 1--20. Available online as latest article  \doi{10.1080/00031305.2023.2186952}
}
\author{ Jan Graffelman (jan.graffelman@upc.edu) }
\seealso{
\code{\link{wAddPCA},\link{ipSymLS},\link{Keller}}
}
\examples{
data(HeartAttack)
X <- HeartAttack[,1:7]
X[,7] <- log(X[,7])
colnames(X)[7] <- "logPR"
R <- cor(X)
W <- matrix(1, 7, 7)
diag(W) <- 0
out.sym <- FitRDeltaQSym(R, W, eps=1e-6) 
Rhat <- out.sym$Rhat
}
\keyword{optimize}
\keyword{multivariate}
