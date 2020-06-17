print.gsda.result <-
function(gsda.result)
{
  gsda.res=unclass(gsda.result)
  gsda.res=as.data.frame(gsda.res)
  rownames(gsda.res)<-NULL
  print(gsda.res[,-2])
}
