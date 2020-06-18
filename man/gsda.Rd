\name{gsda}
\alias{gsda}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Gene-Set Distance Analsysis (GSDA)
}
\description{
This function impletements the gene-set distance analysis (GSDA) of omic data by generalizing distance correlations to evaluate the association of each of a series gene sets with numeric, categorical, and censored event-time variables.
}
\usage{
gsda(omic.data, clin.data, vset.data, clin.vars, omic.dist, clin.dist)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{omic.data}{The genomic data matrix with features as rows and subjects as columns.  The column names of omic.data are assumed to be observation identifiers.  The gsda function calls the function prep.gsda to merge omic.data (by column name) and clin.data (by the column named "ID") before performing the GSDA procedure.}
  \item{clin.data}{A data frame of clinical data. Each row is a subject and each column is a variable.  The "ID" column of clin.data includes observation identifiers. The gsda function calls the function prep.gsda to merge omic.data (by column name) and clin.data (by the column named "ID") before performing the GSDA procedure.  }
  \item{vset.data}{Variable set data.  Each row assigns a variable (column named vID) to a variable set (column named vset).}
  \item{clin.vars}{Column name(s) of clinical variable(s) to be associated with the gene-sets.}
  \item{omic.dist}{The distance metric for omic data, may be "oe" (overall Euclidean), "me" (marginal Euclidean), "om" (overall Manhattan), or "mm" (marginal Manhattan)}
  \item{clin.dist}{The distance metric for clinical data, may be "oe" (overall Euclidean), "me" (marginal Euclidean), "om" (overall Manhattan), or "mm" (marginal Manhattan)}
}
\details{
This function performs the GSDA method described by Cao and Pounds (2020) through generalizing distance correlations to evaluate the association of a gene set with categorical and censored event-time variables. The distance matrices are centered by U-centering and distance correlation is the inner product of the two U-centered distance matrices over the squared of inner product of each of the two U-centered distance matrices. The distance correlation t-statistics asymptotically follows a t-distribution with n*(n-3)/2 degree of freedom according to Zhu et al. (2020). 
}
\value{
A data.frame with the following columns:
\item{vset}{The name of variable set (gene-set).}
\item{vIDs}{The list of variables in the variable set (gene-set).}
\item{dCor}{The distance association statistics for the variable set.}
\item{p.vset}{The p-value.}
\item{comp.time}{Computation time for each set.}
}
\references{
Cao X and Pounds S (2020)  Gene-Set Distance Associations (GSDA): A Powerful Tool for Gene-Set Association Analysis.

Zhu C, Yao S, Zhang X and Shao X (2020) Distance-based and RKHS-based Dependence Metrics in High Dimension. arXiv:1902.03291
}
\author{Xueyuan Cao \email{xcao12@uthsc.edu} and Stanley Pounds \email{stanley.pounds@stjude.org}
}
\seealso{\code{\link{GSDA}}}
\examples{
data(target.aml.clin)
data(target.aml.expr)
data(kegg.ml.gsets)
res=gsda(target.aml.expr,
         target.aml.clin,
         kegg.ml.gsets,
         "Chloroma","oe","ct")
}
\keyword{association;}
\keyword{enrichment;}
\keyword{gene}
