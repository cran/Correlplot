\name{HeartAttack}
\alias{HeartAttack}
\docType{data}
\title{
  Myocardial infarction or Heart attack data
}
\description{
 Data set consisting of 101 observations of patients who suffered a heart attack.   
}
\usage{data("HeartAttack")}
\format{
  A data frame with 101 observations on the following 8 variables.
  \describe{
    \item{\code{Pulse}}{Pulse}
    \item{\code{CI}}{Cardiac index}
    \item{\code{SI}}{Systolic index}
    \item{\code{DBP}}{Diastolic blood pressure}
    \item{\code{PA}}{Pulmonary artery pressure}
    \item{\code{VP}}{Ventricular pressure}
    \item{\code{PR}}{Pulmonary resistance}
    \item{\code{Status}}{Deceased or survived}
  }
}
\source{Table 18.1, (Saporta 1990, pp. 452--454)}
\references{Saporta, G. (1990) Probabilit\'{e}s analyse des donn\'{e}es et statistique. Paris,\'{E}ditions technip}
\examples{
data(HeartAttack)
str(HeartAttack)
}
\keyword{datasets}
