\name{cat.dist}
\alias{cat.dist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Distance of categorical variable
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
This function calculates distance function for a categorical variable, resulting a nxn square matrix with ij^{th} element to be a indicator function of equality of category membership for i^{th} and j^{th} subject. The distance between each pair of individuals in the same category is zero and the distance between any pair of individuals in different categories is 1.  
}
\value{
A result matrix with nrow and ncol equal to the length of X.
}
\references{
Cao X, Sahr N and Pounds S (2020)  Gene-Set Distance Associations (GSDA): A Powerful Tool for Gene-Set Association Analysis.
}
\author{
Xueyuan Cao <xcao12@uthsc.edu> and Stanley.pounds <stanley.pounds@stjude.org>
}
\seealso{
surv.dist
}
\examples{
data(target.aml.clin)
cd=cat.dist(target.aml.clin$Chloroma)
cd[1:5,1:5]
}