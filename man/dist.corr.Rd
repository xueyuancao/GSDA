\name{dist.corr}
\alias{dist.corr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Distance Correlattion
}
\description{
Calculate the distance correlation for a gene set's omic data matrix with another variable.
}
\usage{
dist.corr(X, Y, x.dist = "me", y.dist = "me")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{The omic numeric data matrix with subjects as rows and variables as columns.  Note this is the TRANSPOSE of how some omic data analysis packages represent omic data and how the omic data is represented in the gsda function of this package.}
  \item{Y}{Numeric data matrix, vector, or data.frame.  The rows of X and rows of Y must represent the same set of subjects in the same order.}
  \item{x.dist}{The distance metric for omic data (X), may be "oe" (overall Euclidean), "me" (marginal Euclidean), "om" (overall Manhattan), or "mm" (marginal Manhattan).}
  \item{y.dist}{The distance metric for clinical data, may be "oe" (overall Euclidean), "me" (marginal Euclidean), "om" (overall Manhattan), or "mm" (marginal Manhattan), same options as for X}
}
\details{
The function calculates distance matrix for X and Y using one of the four methods "oe" (overall Euclidean), "me" (marginal Euclidean), "om" (overall Manhattan), or "mm" (marginal Manhattan). Then, the distance matrices are centered by U-centering and distance correlation is calculated as the inner product of the two U-centered distance matrices over the squared of inner product of each of the two U-centered distance matrices. The distance correlation t-statistics follows a t-distribution with n*(n-3)/2 degree of freedom according to Zhu et al.(2020).
}
\value{
A list with the following components:
\item{odCor}{overall distance correlation statistic}
\item{t.odCor}{t-stat for overall distance correlation statistic}
\item{p.odCor}{p-value for overall distance correlation statistic}
\item{dCor}{distance-based correlation matrix for each pair of variables.}
\item{t.dCor}{t-stat for distance-based correlation matrix}
\item{p.dCor}{p-value for distance-based correlation matrix}
\item{X}{echo input data matrice X}
\item{Y}{echo input data matrice Y}
\item{x.dist}{echo input distance metric for X}
\item{y.dist}{echo input distance metric for Y}
}
\references{
Cao X and Pounds S (2020)  Gene-Set Distance Associations (GSDA): A Powerful Tool for Gene-Set Association Analysis.

Zhu C, Yao S, Zhang X and Shao X (2020) Distance-based and RKHS-based Dependence Metrics in High Dimension. arXiv:1902.03291
}
\author{Xueyuan Cao \email{xcao12@uthsc.edu} and Stanley Pounds \email{stanley.pounds@stjude.org}
}
\seealso{\code{\link{best.dist.corr}}}
\examples{
data(target.aml.clin)
data(target.aml.expr)
target.aml.expr=sqrt(target.aml.expr)
target.aml.expr=t(target.aml.expr)

dc.chl=dist.corr(target.aml.expr,
                 target.aml.clin$Chloroma,
                 "oe","ct")
}