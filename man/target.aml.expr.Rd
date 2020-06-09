\name{target.aml.expr}
\alias{target.aml.expr}
\docType{data}
\title{
RNA-seq expression from the AML TARGET project
}
\description{
A matrix of RNA-seq gene expression values for 123 pediatric AML patients from the TARGET project for genes in the KEGG AML and CML pathways
}
\usage{data("target.aml.expr")}
\format{
  The format is:
 num [1:94, 1:123] 9.69 4.54 103.15 7.79 235.79 ...
 - attr(*, "dimnames")=List of 2
  ..$ : chr [1:94] "ENSG00000002330" "ENSG00000007968" "ENSG00000033327" "ENSG00000051382" ...
  ..$ : chr [1:123] "TARGET-20-PABGKN" "TARGET-20-PACDZR" "TARGET-20-PADDXZ" "TARGET-20-PADYIR" ...
}
\details{
A RNA-seq dataset with expression levels in 94 ensemble gene identifiers for 123 pediatric AML patients
}
\source{
\url{ftp://caftpd.nci.nih.gov/pub/OCG-DCC/TARGET/AML/mRNA-seq/L3/expression/}
}
\examples{
data(target.aml.expr)
}
\keyword{datasets}
