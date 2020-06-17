\name{surv.dist}
\alias{surv.dist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Distance of a Survival Endpoint
}
\description{
A function to calculate the distance for a survivial endpoint. 
}
\usage{
surv.dist(stime.evnt)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{stime.evnt}{
A data frame with time in first column and censor in second column.
}
}
\details{
This function calculates the distance matrix for a censored event-time variable. The calculation is based on the formula in Section 2.4 of Cao and Pounds (2020). The distance metric for censored event-time data is based on the rank-based association statistic for this type of data proposed by Jung et al (2005). 
}
\value{
A square matrix with nrow and ncol equal to the nrow of stime.evnt.  Entry (i,j) of the result matrix gives the survival distance between subjects represented in rows i and j of stime.evnt.
}
\references{
Cao X and Pounds S (2020)  Gene-Set Distance Associations (GSDA): A Powerful Tool for Gene-Set Association Analysis. \cr \cr 
Jung SH, Owzar K, and George SL (2005) A mutiple testing procedure to associate gene expression levels with survival. Statistics in Medicine 24: 3077-88.
}
\author{
Xueyuan Cao <xcao12@uthsc.edu> and Stanley Pounds <stanley.pounds@stjude.org>
}

\seealso{\code{\link{cat.dist}}}
\examples{
data(target.aml.clin)
srv.dist=surv.dist(target.aml.clin[,c("efs.time","efs.evnt")])
}

