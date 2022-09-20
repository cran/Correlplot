\name{Kernels}
\alias{Kernels}
\docType{data}
\title{
  Wheat kernel data
}
\description{
  Wheat kernel data set taken from the UCI Machine Learning Repository
}
\usage{data("Kernels")}
\format{
  A data frame with 210 observations on the following 8 variables.
  \describe{
    \item{\code{area}}{Area of the kernel}
    \item{\code{perimeter}}{Perimeter of the kernel}
    \item{\code{compactness}}{Compactness (C = 4*pi*A/P^2)}
    \item{\code{length}}{Length of the kernel}
    \item{\code{width}}{Width of the kernel}
    \item{\code{asymmetry}}{Asymmetry coefficient}
    \item{\code{groove}}{Length of the groove of the kernel}
    \item{\code{variety}}{Variety (1=Kama, 2=Rosa, 3=Canadian)}
  }
}
\source{
    https://archive.ics.uci.edu/ml/datasets/seeds
}
\references{
M. Charytanowicz, J. Niewczas, P. Kulczycki, P.A. Kowalski, S. Lukasik, S. Zak, A Complete Gradient Clustering Algorithm for Features Analysis of X-ray Images. in: Information Technologies in Biomedicine, Ewa Pietka, Jacek Kawa (eds.), Springer-Verlag, Berlin-Heidelberg, 2010, pp. 15-24.
}
\examples{
data(Kernels)
}
\keyword{datasets}
