\name{wAddPCA}
\alias{wAddPCA}
\title{
  Low-rank matrix approximation by weighted alternating least squares
}
\description{
  Function \code{wAddPCA} calculates a weighted least squares approximation of low rank to a given matrix. 
}
\usage{
wAddPCA(x, w = matrix(1, nrow(x), ncol(x)), p = 2, add = "all", bnd = "opt",
        itmaxout = 1000, itmaxin = 1000, epsout = 1e-06, epsin = 1e-06,
	verboseout = TRUE, verbosein = FALSE)
}
\arguments{
  \item{x}{The data matrix to be approximated}
  \item{w}{The weight matrix}
  \item{p}{The dimensionality of the low-rank solution (2 by default)}
  \item{add}{The additive adjustment to be employed. Can be "all" (default), "nul" (no adjustment),
             "one" (adjustment by a single scalar), "row" (adjustment by a row) or "col" (adjustment by a
	     column).}
  \item{bnd}{Can be "opt" (default), "all", "row" or "col".}
  \item{itmaxout}{Maximum number of iterations for the outer loop of the algorithm}
  \item{itmaxin}{Maximum number of iterations for the inner loop of the algorithm}
  \item{epsout}{Numerical criterion for convergence of the outer loop}
  \item{epsin}{Numerical criterion for convergence of the inner loop}
  \item{verboseout}{Be verbose on the outer loop iterations}
  \item{verbosein}{Be verbose on the inner loop iterations}
}
\value{A list object with fields:
\item{a}{The left matrix (A) of the factorization X = AB'}
\item{b}{The right matrix (B) of the factorization X = AB'}
\item{z}{The product AB'}
\item{f}{The final value of the loss function}
\item{u}{Vector for rows used to construct rank 1 weights}
\item{v}{Vector for columns used to construct rank 1 weights}
\item{p}{The vector with row adjustments}
\item{q}{The vector with column adjustments}
\item{itel}{Iterations needed for convergence} 
\item{delta}{The additive adjustment}
\item{y}{The low-rank approximation to \code{x}}
}
\references{
Graffelman, J. and De Leeuw, J. (2023) Improved approximation and visualization of the correlation matrix. The American Statistician pp. 1--20. Available online as latest article  \doi{10.1080/00031305.2023.2186952}

https://jansweb.netlify
}
\author{jan@deleeuwpdx.net}
\seealso{
\code{\link{ipSymLS}}
}
\examples{
data(HeartAttack)
X <- HeartAttack[,1:7]
X[,7] <- log(X[,7])
colnames(X)[7] <- "logPR"
R <- cor(X)
W <- matrix(1, 7, 7)
diag(W) <- 0
Wals.out <- wAddPCA(R, W, add = "nul", verboseout = FALSE) 
Rhat <- Wals.out$y
}
\keyword{optimize}
\keyword{multivariate}
