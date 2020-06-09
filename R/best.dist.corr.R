best.dist.corr <-
function(X,  # numeric data matrix, rows for subjects and columns for variables
                        Y,  # numeric data matrix, vector, or data.frame
                        x.dist="oe", # distance method for X, may be "oe" for overall Euclidean, "me" for marginal Euclidean, "om" for overall Manhattan, "mm" for marginal Manhattan, "ct" for categorical, or "st" for censored survival time
                        y.dist="oe") # distance method for Y, same options as for X
  
{
  dc.all=dist.corr(X,Y,x.dist,y.dist)
  if (!is.matrix(X))
  {
    warning("X is not a matrix.  No search performed.  Returning dist.corr result.")
    return(dc.all)
  }
  
  k=ncol(X)
  if (k ==1)
  {
    warning("X has only one column.  No search performed.  Returning dist.corr result.")
    return(dc.all)
  }
  
  if (!is.element(x.dist,c("oe","me","om","mm")))
  {
    warning("The distance measure for X must be Euclidean or Manhattan to perform a search.  Returning dist.corr result.")
    return(dc.all)
  }
  
  # Now find best subset by backwards elimination
  out.res=cbind(drop.var=0:(k-1),
                nlog10p=NA)
  out.res[1,2]=-log10(dc.all$p.odCor)
  
  rX=X                  # remaining X matrix
  r.indx=(1:k)          # remaining column index
  rk=k                  # number of remaining columns
  for (i in 1:(k-1))    # Loop over variables to drop
  {
    best.p=1.1                        # Initialize the best p-value among remaining variables to drop         
    drop.indx=NA                      # Initialize index of best variable to drop
    for (j in 1:rk)                   # Loop over remaining variables
    {
      tX=rX[,-j]                        # Drop variable j from the remaining matrix
      test.res=dist.corr(tX,Y,          # Compute the distance correlation
                         x.dist,
                         y.dist)
      
      if (test.res$p.odCor<best.p)      # if this result has a smaller p-value
      {                                 
        best.p=test.res$p.odCor           # update the best p-value
        drop.indx=j                       # update the index of the best p-value
      }
    }                                 # End loop over remaining variables
    out.res[i+1,1]=r.indx[drop.indx]  # Record the index of the variable dropped in this round 
    out.res[i+1,2]=-log10(best.p)     # Record the minus log10 of the best p-value
    r.indx=r.indx[-drop.indx]         # Remove the dropped variable from the list of remaining variables
    rX=rX[,-drop.indx]                # Remove the dropped variable from the matrix of remaining variables
    rk=length(r.indx)                 # Update the number of remaining variables
  }                  # End loop over variables to drop
  
  best.indx=which.max(out.res[,2])
  drop.indx=out.res[1:best.indx,1]
  rX=X[,-drop.indx]
  
  #############################################
  # result object: a list with the following components
  res=list(rX=rX,                        # reduced X matrix
           best.res=out.res[best.indx,], # best result by backward elimination
           all.res=out.res,              # all backward elimination results
                                         # the first column has the index of the column of X that was dropped
                                         # the second column has the negative log10 p-value of the resulting X matrix
    
           X=X,              # echoes input X
           Y=Y)              # echoes input Y
  
  class(res)="bdc"
  
  return(res)
}
