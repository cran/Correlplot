\name{aircraft}
\alias{aircraft}
\docType{data}
\title{Characteristics of aircraft}
\description{Four variables registered for 21 types of aircraft.}
\usage{data("aircraft")}
\format{
  A data frame with 21 observations on the following 4 variables.
  \describe{
    \item{\code{SPR}}{specific power}
    \item{\code{RGF}}{flight range factor}
    \item{\code{PLF}}{payload}
    \item{\code{SLF}}{sustained load factor}
  }
}
\source{Gower and Hand, Table 2.1}
\references{Gower, J.C. and Hand, D.J. (1996) \emph{Biplots}, Chapman & Hall, London}
\examples{
data(aircraft)
str(aircraft)
}
\keyword{datasets}
