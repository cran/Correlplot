\name{pfa}
\alias{pfa}
\title{Principal factor analysis}
\description{
Program \code{pfa} performs (iterative) principal factor analysis, which
is based on the computation of eigenvalues of the reduced correlation matrix.
}
\usage{
pfa(X, option = "data", m = 2, initial.communality = "R2", crit = 0.001, verbose = FALSE)
}
\arguments{
  \item{X}{A data matrix or correlation matrix}
  \item{option}{Specifies the type of matrix supplied by argument
    \code{X}. Values for \code{option} are \code{data}, \code{cor} or
    \code{cov}. \code{data} is the default.}
  \item{m}{ The number of factors to extract (2 by default)}
  \item{initial.communality}{ Method for computing initial
    communalites. Possibilities are \code{R2}  or \code{maxcor}.}
  \item{crit}{ The criterion for convergence. The default is
    \code{0.001}. A smaller value will require more iterations before
    convergence is reached.}
  \item{verbose}{ When set to \code{TRUE}, additional numerical output
    is shown.}
}
\value{
  \item{Res }{ Matrix of residuals }
  \item{Psi }{ Diagonal matrix with specific variances }
  \item{La }{ Matrix of loadings }
  \item{Shat }{ Estimated correlation matrix}
  \item{Fs }{ Factor scores} 
}
\references{
  Mardia, K.V., Kent, J.T. and Bibby, J.M. (1979) Multivariate analysis.
  
  Rencher, A.C. (1995) Methods of multivriate analysis.

  Satorra, A. and Neudecker, H. (1998) Least-Squares Approximation of
  off-Diagonal Elements of a Variance Matrix in the Context of Factor
  Analysis. Econometric Theory 14(1) pp. 156--157.
}
\author{ Jan Graffelman (jan.graffelman@upc.edu) }
\seealso{\code{\link{princomp}}}
\examples{
   X <- matrix(rnorm(100),ncol=2)
   out.pfa <- pfa(X)
#  based on a correlation matrix
   R <- cor(X)
   out.pfa <- pfa(R,option="cor")
}
\keyword{multivariate}

