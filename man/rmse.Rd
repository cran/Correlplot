\name{rmse}
\alias{rmse}
\title{
  Calculate the root mean squared error
}
\description{
  Program \code{rmse} calculates the RMSE for matrix using all elements or all off-diagonal matrix.
}
\usage{
rmse(R, Rhat, option = "offd", verbose = FALSE)
}
\arguments{
  \item{R}{The original matrix}
  \item{Rhat}{The approximating matrix}
  \item{option}{Use all elements (\code{option="full"}) or all off-diagonal elements (\code{option="offd"})}
  \item{verbose}{Print output (\code{verbose=TRUE}) or not (\code{verbose=FALSE})}
}
\value{the calculated rmse}
\author{ Jan Graffelman (jan.graffelman@upc.edu) }
\examples{
data(banknotes)
X <- as.matrix(banknotes[,1:6])
R <- cor(X)
out.sd <- eigen(R)
V <- out.sd$vectors
Dl <- diag(out.sd$values)
V2 <- V[,1:2]
D2 <- Dl[1:2,1:2]
Rhat <- V2\%*\%D2\%*\%t(V2)
rmse(R,Rhat,option="offd")
}
\keyword{misc}
