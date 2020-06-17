\name{kegg.ml.gsets}
\alias{kegg.ml.gsets}
\docType{data}
\title{
KEGG gene set data for the AML and CML pathways
}
\description{
A data set with the list of ensemble gene identifiers for the acute myeloid leukemia (AML) and chronic myeloid leukemia pathways as defined in the KEGG pathway database
}
\usage{data("kegg.ml.gsets")}
\format{
  A data frame with 128 rows describing the pairings between the following 2 variables.
  \describe{
    \item{\code{vset}}{KEGG pathway name}
    \item{\code{vID}}{Ensemble gene (ENSG) identifier}
  }
}
\details{A dataset with assignments of ensemble gene identifiers (ENSG) to KEGG pathway names
}
\source{
\url{http://software.broadinstitute.org/gsea/msigdb/cards/KEGG_CHRONIC_MYELOID_LEUKEMIA} \cr \cr
\url{http://software.broadinstitute.org/gsea/msigdb/cards/KEGG_ACUTE_MYELOID_LEUKEMIA} \cr \cr
Merged with information from \url{ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz}
to translate gene symbols to Ensemble gene (ENSG) identifiers
}
\examples{
data(kegg.ml.gsets)
}
\keyword{datasets}
