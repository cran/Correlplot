\name{banknotes}
\alias{banknotes}
\docType{data}
\title{
  Swiss banknote data
}
\description{
  The Swiss banknote data consist of six measures taken on 200 banknotes, of which 100 are counterfeits, and 100 are normal. 
}
\usage{data("banknotes")}
\format{
  A data frame with 200 observations on the following 7 variables.
  \describe{
    \item{\code{Length}}{Banknote length}
    \item{\code{Left}}{Left width}
    \item{\code{Right}}{Right width}
    \item{\code{Bottom}}{Bottom margin}
    \item{\code{Top}}{Top margin}
    \item{\code{Diagonal}}{Length of the diagonal of the image}
    \item{\code{Counterfeit}}{0 = normal, 1 = counterfeit}
  }
}
\references{
Weisberg, S. (2005) Applied Linear Regression. Third edition. John Wiley & Sons, New Jersey.
}
\examples{
data(banknotes)
}
\keyword{datasets}
