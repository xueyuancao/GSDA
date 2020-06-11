## ----loadLibs, include = FALSE-----------------------
#library(MASS)
library(GSDA)
#library(ggplot2)
library(knitr)
data(target.aml.clin)
data(target.aml.expr)
data(kegg.ml.gsets)
opts_chunk$set(
  comment = "#>",
  collapse = TRUE,
  digits = 3,
  tidy = FALSE,
  background = "#FFFF00",
  fig.align = 'center',
  warning = FALSE,
  message = FALSE
  )
options(width = 55, digits = 3)
#theme_set(theme_bw())
getInfo <- function(what = "Suggests") {
  text <- packageDescription("GSDA")[what][[1]]
  text <- gsub("\n", ", ", text, fixed = TRUE)
  text <- gsub(">=", "$\\\\ge$", text, fixed = TRUE)
  eachPkg <- strsplit(text, ", ", fixed = TRUE)[[1]]
  eachPkg <- gsub(",", "", eachPkg, fixed = TRUE)
  #out <- paste("\\\**", eachPkg[order(tolower(eachPkg))], "}", sep = "")
  #paste(out, collapse = ", ")
  length(eachPkg)
}

## ----Chloroma----------------------------------------
chl.res=gsda(sqrt(target.aml.expr),
             target.aml.clin,
             kegg.ml.gsets,
             "Chloroma","oe","ct")
chl.res

## ----logWBC------------------------------------------
wbc.res=gsda(sqrt(target.aml.expr),
             target.aml.clin,
             kegg.ml.gsets,
             "logWBC","oe","oe")
wbc.res


## ----efs---------------------------------------------
efs.res=gsda(sqrt(target.aml.expr),
             target.aml.clin,
             kegg.ml.gsets,
             c("efs.time","efs.evnt"),"oe","st")
efs.res


## ----install, eval = FALSE---------------------------
#  devtools::install_github('xueyuancao/GSDA', dependencies = c("Depends", "Suggests"))

