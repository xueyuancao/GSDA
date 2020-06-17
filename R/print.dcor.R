print.dcor <-
function(dcor.result) # result of dist.corr
{
  cat("--------- Distance Correlation Result -------------\n")
  
  cat(paste("Distance Metric for X:",dist.metric.name(dcor.result$x.dist), "\n"))
  
  if (is.element(dcor.result$x.dist,c("oe","me","om","mm")))
  {
    if (is.vector(dcor.result$X))
      dcor.result$X=matrix(dcor.result$X,ncol=1)
    cat(paste("  Dimension of X: ",nrow(dcor.result$X),
                ",",ncol(dcor.result$X), "\n"))    
  }

  
  if (dcor.result$x.dist=="ct")
    cat(paste("  Number of X observations:",length(dcor.result$X),
                "; Number of X categories:",length(unique(dcor.result$X)), "\n"))
  
  if (dcor.result$x.dist=="st")
    cat(paste("  Number of X observations:",nrow(dcor.result$X),
                "; Number of X events: ",sum(dcor.result$X[,2]>0), "\n"))
  
  cat(paste("Distance Metric for Y:",dist.metric.name(dcor.result$y.dist), "\n"))
  
  if (is.element(dcor.result$y.dist,c("oe","me","om","mm")))
  {
    if (is.vector(dcor.result$Y))
      dcor.result$Y=matrix(dcor.result$Y,ncol=1)
    cat(paste("  Dimension of Y: ",nrow(dcor.result$Y),
                ",",ncol(dcor.result$Y), "\n"))
  }
    
  
  if (dcor.result$y.dist=="ct")
    cat(paste("  Number of Y observations:",length(dcor.result$Y),
                "; Number of Y categories:",length(unique(dcor.result$Y)), "\n"))
  
  if (dcor.result$y.dist=="st")
    cat(paste("  Number of Y observations:",nrow(dcor.result$Y),
                "; Number of Y events: ",sum(dcor.result$Y[,2]>0), "\n"))

  cat(paste("Overall Distance Correlation:",dcor.result$odCor, "\n"))
  cat(paste("Overall Distance Correlation t-statistic:",dcor.result$t.odCor, "\n"))
  cat(paste("Overall Distance Correlation p-value:",dcor.result$p.odCor, "\n"))
}
