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
 Each row contains the expression values of one 94 Ensemble genes for all 123 patients.  Each column contains the expression values of all 94 Ensemble genes for one patient.  The rownames give the Ensemble identifiers for the genes. The columns give the patient identifiers.  
}
\details{
A RNA-seq dataset with expression levels in 94 ensemble gene identifiers for 123 pediatric AML patients
}
\source{
\url{https://target-data.nci.nih.gov/Public/AML/mRNA-seq/L3/expression/}
}
\examples{
data(target.aml.expr)
}
\keyword{datasets}
