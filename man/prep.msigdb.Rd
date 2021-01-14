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
Take a species from MsigDB (https://www.gsea-msigdb.org/gsea/msigdb/index.jsp), extract gene set definiton and prepare a data frame with gene sets and genes to be used as vset.data in the gsda function.
}
\value{
A two-column data.frame with the columns vset and vID
}
\references{
Cao X and Pounds S (2021)  Gene-Set Distance Associations (GSDA): A Powerful Tool for Gene-Set Association Analysis.
}
\author{Xueyuan Cao \email{xcao12@uthsc.edu} and Stanley Pounds \email{stanley.pounds@stjude.org}
}
\examples{
\donttest{
gsets=prep.msigdb()
head(gsets)
}
}
\keyword{Gene}
\keyword{Set}
