\name{prep.msigdb}
\alias{prep.msigdb}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Preparation of MSigDB for GSDA
}
\description{
This function prepares the gene sets of a species in MsigDB for gene-set distance analysis.
}
\usage{
prep.msigdb(species = "Homo sapiens", vset = "gs_name", vID = "gene_symbol")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{species}{Name of species in MSigDB.}
  \item{vset}{Name of MSigDB column to use as vset in gsda, default is "gs_name".}
  \item{vID}{Name of MSigDB column to use as vID in gsda, default is "gene_symbol".}
}
\details{
Take a species from MsigDB, extract gene set definiton and prepare a data frame with gene sets and genes.
}
\value{
A two-column data.frame with the columns vset and vID
}
\references{
Cao X, Sahr N and Pounds S (2020)  Gene-Set Distance Associations (GSDA): A Powerful Tool for Gene-Set Association Analysis.
}
\author{
Xueyuan Cao <xcao12@uthsc.edu> and Stanley.pounds <stanley.pounds@stjude.org>
}
\examples{
gsets=prep.msigdb()
head(gsets)
}
\keyword{Gene}
\keyword{Set}
