cat.dist <-
function(X) # X = vector of category designations
  
{
  n=length(X)                             # number of observations = number of rows of X
  d=matrix(0,n,n)                         # initialize the matrix
  for (i in 1:(n-1))                      # loop over matrix rows
  {
    for (j in (i+1):n)                    # loop over matrix columns
    {
      temp=as.numeric((X[i]!=X[j]))         # compute distance between subjects i and j
      d[i,j]=temp                           # populate matrix elements  
      d[j,i]=temp          
    }                                     # end loop over columns
  }                                       # end loop over rows  
  tbl=table(X)                            # tabulate unique values of X
  uniq.x=names(tbl)[rev(order(tbl))]      # order unique names by decreasing prevalence 
  attr(d,"vnames")=paste(uniq.x,          # create the name for the result object
                         collapse="|")
  return(d)                               # return the result
}
