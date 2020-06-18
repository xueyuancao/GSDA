\name{uc.dist}
\alias{uc.dist}
\title{
U-centered Distance Matrix
}
\description{
The function calculates the U-centered distance matrix for a variable.
}
\usage{
uc.dist(X, dmeth = "me")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{vector, matrix, or data.frame to compute a distance matrix}
  \item{dmeth}{Distance method to use, options include "oe" for overall Euclidean, "me" for marginal Euclidean, "om" for overall Manhattan, "mm" for marginal Manhattan, "ct" for categorical, and "st" for censored survival time.}
}
\details{
A distance matrix is first calculated for a scale or data frame of a variable. The distance matrix is then centered according to U-centering formula on page 6 of arXiv 1902.03291 paper.
}
\value{
For distance methods "oe", "om", "ct", and "st", one matrix of overall distances computed using data from all variables.  For distance methods "me" and "mm", an array of distance matrices, one distance matrix per variable.
}
\references{
Cao X and Pounds S (2020)  Gene-Set Distance Associations (GSDA): A Powerful Tool for Gene-Set Association Analysis.

Zhu C, Yao S, Zhang X and Shao X. Distance-based and RKHS-based Dependence Metrics in High Dimension. arXiv:1902.03291
}
\author{Xueyuan Cao \email{xcao12@uthsc.edu} and Stanley Pounds \email{stanley.pounds@stjude.org}
}

\seealso{\code{\link{U.center}},  \code{\link{prep.gsda}}}
\examples{
data(target.aml.expr)
target.aml.expr=sqrt(target.aml.expr)
target.aml.expr=t(target.aml.expr)
oe.dist=uc.dist(target.aml.expr,"oe") # overall Euclidean
}
