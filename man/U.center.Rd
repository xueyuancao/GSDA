\name{U.center}
\alias{U.center}
\title{U Centering}
\description{
Center the distance matrix to an comparable scale.
}
\usage{
U.center(d)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d}{A numeric square data matrix}
}
\details{
This funtion Centers the distance matrix according to U-centering formula on page 6 of arXiv 1902.03291 paper
}
\value{A centered data matrix}
\references{
Zhu C, Yao S, Zhang X and Shao X. Distance-based and RKHS-based Dependence Metrics in High Dimension. arXiv:1902.03291
}
\author{
Cao X, Sahr N and Pounds S (2020)  Gene-Set Distance Associations (GSDA): A Powerful Tool for Gene-Set Association Analysis.
}

\seealso{
uc.dist
}
\examples{
data(target.aml.clin)
cd=cat.dist(target.aml.clin$Chloroma)
ud=U.center(cd)
ud[1:5,1:5]
}
