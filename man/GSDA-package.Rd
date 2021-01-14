\name{GSDA-package}
\alias{GSDA-package}
\alias{GSDA}
\docType{package}
\title{
\packageTitle{GSDA}
}
\description{
\packageDescription{GSDA}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{GSDA}
\packageIndices{GSDA}
}
\author{
\packageAuthor{GSDA}

Maintainer: \packageMaintainer{GSDA}
}
\references{
Cao X and Pounds S (2021)  Gene-Set Distance Associations (GSDA): A Powerful Tool for Gene-Set Association Analysis.
}
\keyword{ package }
\examples{
data(target.aml.clin)
data(target.aml.expr)
data(kegg.ml.gsets)
res=gsda(target.aml.expr,
         target.aml.clin,
         kegg.ml.gsets,
         "Chloroma","oe","ct")
}
