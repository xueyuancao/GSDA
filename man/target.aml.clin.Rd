\name{target.aml.clin}
\alias{target.aml.clin}
\docType{data}
\title{Clinical outcomes for AML TARGET Project}
\description{A dataset with subject identifier, survival time, and death indicator for 123 pediatric AML patients}
\usage{data("target.aml.clin")}
\format{
  A data frame with 123 observations of the following 5 variables.
  \describe{
    \item{\code{ID}}{subject identifier, a character vector}
    \item{\code{Chloroma}}{a character vector}
    \item{\code{logWBC}}{a numeric vector}
    \item{\code{efs.time}}{event-free survival time, a numeric vector}
    \item{\code{efs.evnt}}{event indicator (0 = censored, 1 = event) for efs.time, a numeric vector}
  }
}
\details{A dataset with clinical data for each of 123 pediatric AML patients
}
\source{
obtained from \url{ftp://caftpd.nci.nih.gov/pub/OCG-DCC/TARGET/AML/clinical/harmonized/}
}
\examples{
data(target.aml.clin)
}
\keyword{datasets}
