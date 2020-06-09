uc.dist <-
function(X,          # vector, matrix, or data.frame to compute a distance matrix                     
                 dmeth="me") # distance method to use, options include "oe" for overall Euclidean,
                             # "me" for marginal Euclidean, "om" for overall Manhattan,
                             # "mm" for marginal Manhattan, "ct" for categorical,
                             # and "st" for censored survival time
{

  # List of supported distance method
  dmethods=c("oe",    # overall Euclidean
             "om",    # overall Manhattan 
             "me",    # marginal Euclidean
             "mm",    # margian Manhattan
             "ct",    # categorical distance    
             "st")    # survival time distance

  
  # Identify the selected distance method
  dmtch=charmatch(tolower(dmeth),
                  dmethods)
  
  if (dmtch==0)
    stop("Unsupported distance method (dmeth) specified.  Please specify 'euclidean', 'manhattan', 'group', or 'survival'.")
  
  dmeth=dmethods[dmtch]
  
  
  mc=as.list(match.call())
  
  # Overall Euclidean or Manhattan Distance
  if (is.element(dmtch,1:2))
  {
    if (is.vector(X))             # convert X to a matrix if necessary
      X=matrix(X,length(X),1)
    n=nrow(X)                     # number of observations (subjects or samples) in X matrix
    
    dist.meth="euclidean"         # use Euclidean or Manhattan distance
    if (dmtch==2)                 # as specified by the user
      dist.meth="manhattan"
    
    d=dist(X,method=dist.meth)    # then compute Manhattan distance
    
    d=as.matrix(d)                # convert to a matrix
    d=U.center(d)                 # U-center the distance matrix   
    A=array(NA,dim=c(n,n,1))      # Initialize array representation for later calculations
    A[,,1]=d                      # Represent result in array representation
    attr(A,"vnames")=mc$X         # Assign variable name to result
    return(A)                     # return the result
  }
  
  
  # Categorical distance matrix
  if (dmtch==5)
  {
    n=length(X)
    d=cat.dist(X)                             # compute the categorical distance
    d=d/median(abs(d[d>0]))                   # scale for numerical stability
    d=U.center(d)                             # U-center the categorical distance
    A=array(NA,dim=c(n,n,1))                  # Initialize array format for result
    A[,,1]=d                                  # Represent matrix as an array for later calculations
    attr(A,"vnames")=attr(d,"vnames")         # Assign the variable name attribute the result
    return(A)                                 # return the result
  }
  
  # Survival distance matrix
  if (dmtch==6)
  {
    if (is.null(colnames(X)))                 # create column names for X if none were provided
      colnames(X)=paste0("clm",1:ncol(X))
    n=nrow(X)                                 # number of observations = numbers of rows of X
    A=array(NA,dim=c(n,n,1))                  # an array of distance matrices, one distance matrix for each X variable
    d=surv.dist(X)                            # compute survival distance
    d=d/median(abs(d[d>0]))                   # scale the survival distance for numerical stability
    d=as.matrix(d)                            # convert distance result to a matrix
    d=U.center(d)                             # U-center the survival distance
    A[,,1]=d                                  # distance matrix for survival variable
    attr(A,"vnames")=colnames(X)[1]           # assign column name 1 (survival time variable) as variable name for distance matrix array
    return(A)                                 # return the distance matrix array
  }
  
  # Marginal Euclidean or Manhattan Distance
  if (is.element(dmtch,3:4))
  {
    # Distance Matrix for Numeric Values
    if (is.vector(X))                     # convert X to a matrix if necessary
      X=matrix(X,length(X),1)
    n=nrow(X)                             # number of observations (subjects or samples) in X matrix
    k=ncol(X)                             # number of variables in X matrix
    
    if (is.null(colnames(X)))
      colnames(X)=paste0("clm",1:ncol(X))
    
    dist.meth="euclidean"
    if (dmtch==4)
      dist.meth="manhattan"
    
    A=array(NA,dim=c(n,n,k))              # an array of distance matrices, one distance matrix for each X variable
    for (i in 1:k)                        # Loop over variables
    {
      x=X[,i]                                   # extract variable i
      d=dist(x,method=dist.meth)                # compute distances for variable i
      d=d/median(abs(d[d>0]))                   # scale the distance matrix
      d=as.matrix(d)                            # convert distance result to a matrix
      d=U.center(d)
      A[,,i]=d                                  # assign to distance matrix i of the distance array
    }                                    # End loop over variables
    attr(A,"vnames")=colnames(X)
    
    return(A)                            # Return an array of distance matrices, one distance matrix per X variable
    
  }
  
  
}
