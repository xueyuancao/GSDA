\name{print.gsda.result}
\alias{print.gsda.result}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Print Method for GSDA Result
}
\description{
Print the of result of gene-set distance analysis (gsda function).
}
\usage{
\method{print}{gsda.result}(x,...)
%print.gsda.result(gsda.result)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{result of gene-set distance analysis (gsda function)}
  \item{...}{further arguments passed to or from other methods}
}
\details{
Print the result of gene-set distance analysis to stdout.
}
\value{No return value, called for side effects} 
\references{
Cao X and Pounds S (2021)  Gene-Set Distance Associations (GSDA): A Powerful Tool for Gene-Set Association Analysis.
}
\author{Xueyuan Cao \email{xcao12@uthsc.edu} and Stanley Pounds \email{stanley.pounds@stjude.org}
}
\seealso{\code{\link{print.dcor}}, \code{\link{print.bdc}}}
