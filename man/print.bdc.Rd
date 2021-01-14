\name{print.bdc}
\alias{print.bdc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Print Method for Best Distance Correlation
}
\description{
Print the result of the best distance correlation (best.dist.corr)
}
\usage{
\method{print}{bdc}(x,...)
%print.bdc(bdc.result)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a class of bdc}
  \item{...}{further arguments passed to or from other methods}
}
\details{
Print the summary of result of best distance correlation to stdout.
}
\value{No return value, called for side effects} 
\references{
Cao X and Pounds S (2021)  Gene-Set Distance Associations (GSDA): A Powerful Tool for Gene-Set Association Analysis.
}
\author{Xueyuan Cao \email{xcao12@uthsc.edu} and Stanley Pounds \email{stanley.pounds@stjude.org}
}
\seealso{\code{\link{print.dcor}}, \code{\link{print.gsda.result}}}
