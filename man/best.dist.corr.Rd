\name{best.dist.corr}
\alias{best.dist.corr}
\title{Best Distance Correlattion}
\description{
By sequentially dropping varaibles and recaculating the correlation p values, the function find the best set of omic variable in a gene set with best distance correlation p value.  
}
\usage{
best.dist.corr(X, Y, x.dist = "oe", y.dist = "oe")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{The omic numeric data matrix. Each row is a features. Each column is a subject.}
  \item{Y}{Numeric data matrix, vector, or data.frame.}
  \item{x.dist}{The distance metric for omic data (X), may be "oe" (overall Euclidean), "me" (marginal Euclidean), "om" (overall Manhattan), or "mm" (marginal Manhattan).}
  \item{y.dist}{The distance metric for clinical data, may be "oe" (overall Euclidean), "me" (marginal Euclidean), "om" (overall Manhattan), or "mm" (marginal Manhattan), same options as for X}
}
\details{
With full gene set, each gene is recursively dropped to re-compute the distance correlation till best p value is reached. This sub set reflects the driving genes in the distance correlation. 
}
\value{
A list with the following components:
\item{rX}{reduced X matrix}
\item{best.res}{best result by backward elimination}
\item{all.res}{all backward elimination results: the first column has the index of the column of X that was dropped; the second column has the negative log10 p-value of the resulting X matrix}
\item{X}{echoes input X}
\item{Y}{echoes input Y}
}
\references{
Cao X, Sahr N and Pounds S (2020)  Gene-Set Distance Associations (GSDA): A Powerful Tool for Gene-Set Association Analysis.
}
\author{
Xueyuan Cao <xcao12@uthsc.edu> and Stanley.pounds <stanley.pounds@stjude.org>
}
\seealso{dist.corr}
\examples{
data(target.aml.clin)
data(target.aml.expr)
target.aml.expr=sqrt(target.aml.expr)
target.aml.expr=t(target.aml.expr)

bdc.chl=best.dist.corr(target.aml.expr,
                       target.aml.clin$Chloroma,
                       "oe","ct")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
