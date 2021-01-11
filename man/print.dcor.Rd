\name{print.dcor}
\alias{print.dcor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Print Method for Distance Correlation
}
\description{
Print the summary of result of distance correlation (dist.corr function).
}
\usage{
\method{print}{dcor}(x,...)
%print.dcor(dcor.result)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{result of dist.corr, class dcor}
  \item{...}{further arguments passed to or from other methods}
}
\details{
Print the summary of result of distance correlation to stdout.
}
\references{
Cao X and Pounds S (2020)  Gene-Set Distance Associations (GSDA): A Powerful Tool for Gene-Set Association Analysis.
}
\author{Xueyuan Cao \email{xcao12@uthsc.edu} and Stanley Pounds \email{stanley.pounds@stjude.org}
}
\seealso{\code{\link{print.bdc}}, \code{\link{print.gsda.result}}}
