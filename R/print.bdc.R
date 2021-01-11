print.bdc <-
function(x,...) # result of best.dist.corr
  
{
   cat(paste0("Dimension of original X: ",nrow(x$X),", ",ncol(x$X), "\n"))
   cat(paste0("Dimension of reduced X: ",nrow(x$rX),", ",ncol(x$rX), "\n"))
   cat(paste0("Original minus log 10 p-value: ",x$all.res[1,2], "\n"))
   cat(paste0("Minus log 10 Pseudo p-value for reduced X: ",x$best.res[2], "\n"))
   #NextMethod("print")  
}
