\name{Keller}
\alias{Keller}
\title{
  Program \code{Keller} calculates a rank p approximation to a correlation matrix according to Keller's method.
}
\description{
  Keller's method is based on iterated eigenvalue decompositions that are used to adjust the diagonal of the correlation matrix.
}
\usage{
Keller(R, eps = 1e-06, nd = 2, itmax = 10)
}
\arguments{
  \item{R}{A correlation matrix}
  \item{eps}{Numerical criterion for convergence (default \code{eps=1e-06)}}
  \item{nd}{Number of dimensions used in the spectral decomposition (default \code{nd=2})}
  \item{itmax}{The maximum number of iterations}
}
\value{A matrix containing the approximation to the correlation matrix-}
\references{Keller, J.B. (1962) Factorization of Matrices by Least-Squares. Biometrika, 49(1 and 2) pp. 239--242.

Graffelman, J. and De Leeuw, J. (2023) Improved approximation and visualization of the correlation matrix. The American Statistician pp. 1--20. Available online as latest article  \doi{10.1080/00031305.2023.2186952}
}
\author{ Jan Graffelman (jan.graffelman@upc.edu) }
\seealso{
   \code{\link{ipSymLS}}
}
\examples{
data(Kernels)
R <- cor(Kernels)
Rhat <- Keller(R)
}
\keyword{array}
