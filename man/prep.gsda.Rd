\name{prep.gsda}
\alias{prep.gsda}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
GSDA Data Preparation
}
\description{
A function to prepare omic data, clinical data and variable set into an ordered matched format for GSDA analysis.
}
\usage{
prep.gsda(data.mtx, clin.data, vset.data = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data.mtx}{Numeric data matrix with column names giving subject identifiers.}
  \item{clin.data}{Data.frame with column named "ID" with subject identifiers matching column names of data.mtx.}
  \item{vset.data}{data.frame of variable-set assignments with columns named "vID" for variable identifier and "vset" for name or identifier of a variable set (gene-set).}
}
\details{
The gsda function uses prep.gsda to prepare the omic data matrix, clinical dataframe and variable set (gene set) into ordered and matched format, which is then used for GSDA analysis.
}
\value{
A list with the following components:
\item{omic.data}{data matrix with columns in the same order as clin.data$ID.}
\item{clin.data}{data.frame with ID column in same order as columns of omic.data.}
\item{vset.data}{variable set ordered by name of variable set.}
\item{vset.index}{simple data.frame showing first and last row of vset.data for each variable set.}
}
\references{
Cao X and Pounds S (2020)  Gene-Set Distance Associations (GSDA): A Powerful Tool for Gene-Set Association Analysis.
}
\author{Xueyuan Cao \email{xcao12@uthsc.edu} and Stanley Pounds \email{stanley.pounds@stjude.org}
}
\seealso{\code{\link{gsda}}}
\examples{
data(target.aml.clin)
data(target.aml.expr)
data(kegg.ml.gsets)
gsdaprep=prep.gsda(target.aml.expr,
         target.aml.clin,
         kegg.ml.gsets)
}


