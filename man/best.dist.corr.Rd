\name{best.dist.corr}
\alias{best.dist.corr}
\title{Best Distance Correlattion}
\description{
Use a backward elimination procedure to identify a subset of variables in X most strongly associated with Y according to the distance correlation p-value.  
}
\usage{
best.dist.corr(X, Y, x.dist = "oe", y.dist = "oe")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{The omic numeric data matrix with subjects in rows and variables in columns.  Note that this is the TRANSPOSE of the omic data matrix for some other omic data analysis packages and for the gsda function of this package.}
  \item{Y}{Numeric data matrix, vector, or data.frame with each row representing a subject.  The function assumes the same set of subjects are represent in the same order in X and Y.}
  \item{x.dist}{The distance metric for omic data (X), may be "oe" (overall Euclidean), "me" (marginal Euclidean), "om" (overall Manhattan), or "mm" (marginal Manhattan).}
  \item{y.dist}{The distance metric for clinical data, may be "oe" (overall Euclidean), "me" (marginal Euclidean), "om" (overall Manhattan), or "mm" (marginal Manhattan), same options as for X}
}
\details{
This function computes dist.corr for X and Y.  It then determines which column of X may be dropped to give the smallest p-value in dist.corr.  This process is repeated until X has been reduced to only one variable.  In this way, a dist.corr p-value is obtained after dropping each X variable.  The subset of X variables giving the smallest p-value in this series of analyses is returned with additional result details.}
\value{
A list with the following components:
\item{rX}{reduced X matrix}
\item{best.res}{best result by backward elimination}
\item{all.res}{all backward elimination results: the first column has the index of the column of X that was dropped; the second column has the negative log10 p-value of the resulting X matrix}
\item{X}{echoes input X}
\item{Y}{echoes input Y}
}
\references{
Cao X and Pounds S (2020)  Gene-Set Distance Associations (GSDA): A Powerful Tool for Gene-Set Association Analysis.
}
\author{
Xueyuan Cao <xcao12@uthsc.edu> and Stanley.pounds <stanley.pounds@stjude.org>
}
\seealso{\code{\link{dist.corr}}}
\examples{
data(target.aml.clin)
data(target.aml.expr)
target.aml.expr=sqrt(target.aml.expr)
target.aml.expr=t(target.aml.expr)

bdc.chl=best.dist.corr(target.aml.expr,
                       target.aml.clin$Chloroma,
                       "oe","ct")
}
