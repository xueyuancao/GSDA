prep.gsda <-
function(data.mtx,        # numeric data matrix with column names
                   clin.data,       # data.frame with column named "ID" with identifiers matching column names of data.mtx
                   vset.data=NULL)  # data.frame of variable-set assignments with columns named "vID" for variable identifier and "vset" for name or identifier of a variable set (gene-set)
  
{
  clin.clms=colnames(clin.data)
  if (!any(is.element(clin.clms,"ID")))
    stop("clin.data must be a data.frame with a column named ID with column names of data.mtx.")
  
  # Merge data sets, match order of subjects in columns of data.mtx and rows of clin.data
  id.mtx=colnames(data.mtx)                # ids from data.mtx
  id.grp=clin.data$ID                      # ids from clin.data
  id.mtch=match.order.ids(id.mtx,id.grp)   # function to match order of subjects
  data.mtx=data.mtx[,id.mtch$ind1]         # put columns of data.mtx in that order
  clin.data=clin.data[id.mtch$ind2,]       # put rows of clin.data in that order
  
  # check merging before continuing; stop if something didn't work
  id.mtx=colnames(data.mtx)   # ids from data.mtx
  id.grp=clin.data$ID         # ids from clin.data
  if (any(id.mtx!=id.grp))    # stop if something didn't work
    stop("Error in merging clin.data and data.mtx.")
  
  n=nrow(clin.data) # Now get the sample size

  # Order and index variable set data for faster computing
  vset.index=NULL
  if (!is.null(vset.data))
  {
    vset.data0=vset.data                        # keep original variable set data to double-check later
    ord.dset=order.index.dset(vset.data,"vset") # function to order gene-set assignments by gene-set
    vset.data=ord.dset$dset                     # ordered gene-set data
    vset.index=ord.dset$indx                    # indexes of ordered gene-set data
    nset=nrow(vset.index)                       # number of gene-sets
    
    # Double-check the re-indexed variable set data
    ord0=order(vset.data0$vset,vset.data0$vID)  # order original and re-indexed vset data in the same way
    ord1=order(vset.data$vset,vset.data$vID)
    if (any(vset.data0$vID[ord0]!=vset.data$vID[ord1])||
        any(vset.data0$vset[ord0]!=vset.data$vset[ord1]))
    {
      stop("Error in indexing vset.data.")
    }
    
    # Check that each indexed variable set has only one value for vset variable
    for (i in 1:nset)
    {
      indx=vset.index[i,1]:vset.index[i,2]
      uniq.vset=unique(vset.data$vset[indx])
      if (length(uniq.vset)!=1)
        stop("Error in index vset.data")
    } 
  }
  
  ##############################
  # result: a list with the following components
  res=list(omic.data=data.mtx,    # data matrix with columns in the same order as clin.data$ID
           clin.data=clin.data,   # data.frame with ID column in same order as columns of omic.data
           vset.data=vset.data,   # variable set ordered by name of variable set
           vset.index=vset.index) # simple data.frame showing first and last row of vset.data for each variable set
  
  return(res)
  
}
