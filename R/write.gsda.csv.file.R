write.gsda.csv.file <-
function(gsda.result,out.file)
{
  if (substring(out.file,nchar(out.file)-3)!=".csv")
    stop("out.file must be a .csv file!")
  gsda.res=unclass(gsda.result)
  gsda.res=as.data.frame(gsda.res)
  write.csv(gsda.res,out.file,
            quote=T,row.names=F)
}
