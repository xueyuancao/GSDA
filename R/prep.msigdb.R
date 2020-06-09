prep.msigdb <-
function(species="Homo sapiens", # Name of species in MSigDB
                     vset="gs_name",         # Name of MSigDB column to use as vset in gsda
                     vID="gene_symbol")      # Name of MSigDB column to use as vID in gsda
{
  species.list=msigdbr_show_species()             
  if (!is.element(tolower(species),
                  tolower(species.list)))
  {
    warning("Invalid species request.")
    message("Choose from one of the following:")
    print(species.list)
    stop("Invalid species request.")
  }
  
  res=msigdbr(species=species)
  res=as.data.frame(res)
  if (any(!is.element(c(vset,vID),colnames(res))))
  {
    warning("Invalid specification of vset or vID.")
    print("Valid choices for vset or vID are:")
    print(colnames(res))
    stop("Invalid specification of vset or vID.")
  }
    
  # Returns a two-column data.frame with the columns vset and vID
  res=res[,c("gs_name","gene_symbol")]
  names(res)=c("vset","vID")
  return(res)
}
