print.bdc <-
function(bdc.result,...) # result of best.dist.corr
  
{
   cat(paste0("Dimension of original X: ",nrow(bdc.result$X),", ",ncol(bdc.result$X), "\n"))
   cat(paste0("Dimension of reduced X: ",nrow(bdc.result$rX),", ",ncol(bdc.result$rX), "\n"))
   cat(paste0("Original minus log 10 p-value: ",bdc.result$all.res[1,2], "\n"))
   cat(paste0("Minus log 10 Pseudo p-value for reduced X: ",bdc.result$best.res[2], "\n"))  
}
