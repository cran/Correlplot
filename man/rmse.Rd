\name{rmse}
\alias{rmse}
\title{
  Calculate the root mean squared error
}
\description{
  Program \code{rmse} calculates the RMSE for a matrix approximation.
}
\usage{
rmse(R, Rhat, W = matrix(1, nrow(R), ncol(R)) - diag(nrow(R)),
     verbose = FALSE, per.variable = FALSE)
}
\arguments{
  \item{R}{The original matrix}
  \item{Rhat}{The approximating matrix}
  \item{W}{A symmetric matrix of weights}
  \item{verbose}{Print output (\code{verbose=TRUE}) or not (\code{verbose=FALSE})}
  \item{per.variable}{Calculate the RMSE for the whole matrix (\code{per.variable=FALSE}) or
        for each variable seperately (\code{per.variable=TRUE})}
}
\details{
By default, function \code{rmse} assumes a symmetric correlation matrix as input, together with its approximation. The approximation does not need to be symmetric.
Weight matrix W has to be symmetric. By default, the diagonal is excluded from RMSE calcuations (W = J - I). To include it, specify W = J, that is set W = matrix(1, nrow(R), ncol(R))
}
\value{the calculated rmse}
\references{Graffelman, J. and De Leeuw, J. (2023) Improved approximation and visualization of the correlation matrix. The American Statistician pp. 1--20. \doi{10.1080/00031305.2023.2186952}}
\author{ Jan Graffelman (jan.graffelman@upc.edu) }
\examples{
data(banknotes)
X <- as.matrix(banknotes[,1:6])
p <- ncol(X)
J <- matrix(1,p,p)
R <- cor(X)
out.sd <- eigen(R)
V <- out.sd$vectors
Dl <- diag(out.sd$values)
V2 <- V[,1:2]
D2 <- Dl[1:2,1:2]
Rhat <- V2\%*\%D2\%*\%t(V2)
rmse(R,Rhat,W=J)
}
\keyword{misc}
