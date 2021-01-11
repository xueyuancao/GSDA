print.gsda.result <-
function(x,...)
{
  gsda.res=unclass(x)
  gsda.res=as.data.frame(gsda.res)
  rownames(gsda.res)<-NULL
  print(gsda.res[,-2],...)
}
