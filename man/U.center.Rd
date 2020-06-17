\name{U.center}
\alias{U.center}
\title{U Centering}
\description{
U-center the distance matrix in preparation of computing distance correlations.
}
\usage{
U.center(d)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d}{A square numeric data matrix}
}
\details{
This funtion centers the distance matrix according to U-centering formula on page 6 of arXiv 1902.03291 paper
}
\value{A centered data matrix}
\references{
Zhu C, Yao S, Zhang X and Shao X. Distance-based and RKHS-based Dependence Metrics in High Dimension. arXiv:1902.03291
}
\author{
Cao X and Pounds S (2020)  Gene-Set Distance Associations (GSDA): A Powerful Tool for Gene-Set Association Analysis.
}

\seealso{\code{\link{uc.dist}}}
\examples{
data(target.aml.clin)
cd=cat.dist(target.aml.clin$Chloroma)
ud=U.center(cd)
ud[1:5,1:5]
}
