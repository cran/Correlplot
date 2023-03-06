\name{ipSymLS}
\alias{ipSymLS}
\title{
  Function for obtaining a weighted least squares low-rank approximation of a symmetric matrix
}
\description{
  Function \code{ipSymLS} implements an alternating least squares algorithm that uses both decomposition and block relaxation
  to find the optimal positive semidefinite approxation of given rank p to a known symmetric matrix of order n.
}
\usage{
ipSymLS(target, w = matrix(1, dim(target)[1], dim(target)[2]), ndim = 2,
        init = FALSE, itmax = 100, eps = 1e-06, verbose = FALSE)
}
\arguments{
  \item{target}{Symmetric matrix to be approximated}
  \item{w}{Matrix of weights}
  \item{ndim}{Number of dimensions extracted (2 by default)}
  \item{init}{Initial value for the solution (optional; if supplied should be a matrix of dimensions \code{nrow(target)} by \code{ndim})}
  \item{itmax}{Maximum number of iterations}
  \item{eps}{Tolerance criterion for convergence}
  \item{verbose}{Show the iteration history (\code{verbose=TRUE}) or not (\code{verbose=FALSE})}
}
\value{A matrix with the coordinates for the variables}
\references{De Leeuw, J. (2006) A decomposition method for weighted least squares low-rank approximation of symmetric matrices. Department of Statistics, UCLA. Retrieved from https://escholarship.org/uc/item/1wh197mh

Graffelman, J. and De Leeuw, J. (2023) Improved approximation and visualization of the correlation matrix. The American Statistician pp. 1--20. Available online as latest article  \doi{10.1080/00031305.2023.2186952}
}
\author{deleeuw@stat.ucla.edu}
\examples{
data(banknotes)
R <- cor(banknotes)
W <- matrix(1,nrow(R),nrow(R))
diag(W) <- 0
Fp.als <- ipSymLS(R,w=W,verbose=TRUE,eps=1e-15)
Rhat.als <- Fp.als\%*\%t(Fp.als)
}
\keyword{optimize}
