print.dcor <-
function(x, ...) # result of dist.corr
{
  cat("--------- Distance Correlation Result -------------\n")
  
  cat(paste("Distance Metric for X:",dist.metric.name(x$x.dist), "\n"))
  
  if (is.element(x$x.dist,c("oe","me","om","mm")))
  {
    if (is.vector(x$X))
      x$X=matrix(x$X,ncol=1)
    cat(paste("  Dimension of X: ",nrow(x$X),
                ",",ncol(x$X), "\n"))    
  }

  
  if (x$x.dist=="ct")
    cat(paste("  Number of X observations:",length(x$X),
                "; Number of X categories:",length(unique(x$X)), "\n"))
  
  if (x$x.dist=="st")
    cat(paste("  Number of X observations:",nrow(x$X),
                "; Number of X events: ",sum(x$X[,2]>0), "\n"))
  
  cat(paste("Distance Metric for Y:",dist.metric.name(x$y.dist), "\n"))
  
  if (is.element(x$y.dist,c("oe","me","om","mm")))
  {
    if (is.vector(x$Y))
      x$Y=matrix(x$Y,ncol=1)
    cat(paste("  Dimension of Y: ",nrow(x$Y),
                ",",ncol(x$Y), "\n"))
  }
    
  
  if (x$y.dist=="ct")
    cat(paste("  Number of Y observations:",length(x$Y),
                "; Number of Y categories:",length(unique(x$Y)), "\n"))
  
  if (x$y.dist=="st")
    cat(paste("  Number of Y observations:",nrow(x$Y),
                "; Number of Y events: ",sum(x$Y[,2]>0), "\n"))

  cat(paste("Overall Distance Correlation:",x$odCor, "\n"))
  cat(paste("Overall Distance Correlation t-statistic:",x$t.odCor, "\n"))
  cat(paste("Overall Distance Correlation p-value:",x$p.odCor, "\n"))
  #NextMethod("print")
}
