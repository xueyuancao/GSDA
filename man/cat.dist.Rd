\name{cat.dist}
\alias{cat.dist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Distance for a Categorical Variable
}
\description{
A function to calculate the distance for a categorical variable.
}
\usage{
cat.dist(X)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{vector of category designations}
}
\details{
This function calculates distance function for a categorical variable.  The result is a square n by n matrix in which entry (i,j) has value 1 if entry i and entry j of the input vector X are not equal and entry (i,j) of the result matrix has value 0 if entry i and entry j of the input vector are equal.  The distance between subject i and subject j is zero if the two subjects have the same categorical designation.  The distance between subject i and subject j is one if the two subjects do not have the same categorical designation.
}
\value{
A square matrix with each dimension equal to the length of X.
}
\references{
Cao X and Pounds S (2020)  Gene-Set Distance Associations (GSDA): A Powerful Tool for Gene-Set Association Analysis.
}
\author{
Xueyuan Cao <xcao12@uthsc.edu> and Stanley.pounds <stanley.pounds@stjude.org>
}
\seealso{\code{\link{surv.dist }}}
\examples{
data(target.aml.clin)
cd=cat.dist(target.aml.clin$Chloroma)
cd[1:5,1:5]
}