dist.corr <-
function(X,              # data set 1, maybe a multivariate data matrix, a group label vector, or 2-column survival data matrix (rows for subjects)
                   Y,              # data set 2 may be a multivariate data matrix, a group label vector, or 2-column survival data matrix (rows for subjects)
                   x.dist="me",    # distance metric for X, may be overall Euclidean "oe", overall Manhattan "om", marginal Euclidean "me", categorical "ct", or survival time "st"
                   y.dist="me")    # distance metric for Y, same options as above

{
  ##############################################
  # Check that X and Y have an equal number of observatons
  if (is.vector(X)) nX=length(X)
  else nX=nrow(X)
  
  if (is.vector(Y)) nY=length(Y)
  else nY=nrow(Y)
  
  if (nX!=nY)
    stop("X and Y do not have the same number of observations.")

  
  ###################################################
  # compute U-centered distances for each variable (column) of each matrix
  # \\tilde{a}_{st} defined on page 6 of Zhu et al
  uc.X=uc.dist(X,x.dist)
  uc.Y=uc.dist(Y,y.dist)
  
  ##############################################
  # Sample size and dimension information
  n=nX             # sample size
  kX=dim(uc.X)[3]  # number of X variables
  kY=dim(uc.Y)[3]  # number of Y variables   

  
  ########################################################
  # Compute distance covariances for each pair of variables

  dCov=matrix(NA,kX+kY,kX+kY)  # Initialize the matrix
  colnames(dCov)=c(attr(uc.X,"vnames"),
                   attr(uc.Y,"vnames"))
  rownames(dCov)=colnames(dCov)
  
  # dCov(X,X)
  for (i in 1:kX)  # loop over first dimension of X for pairs of dimensions of X
  {
    for (j in 1:i) # loop over second dimension of X for pairs of dimensions of X
    {
      prd=uc.X[,,i]*uc.X[,,j]                          # compute elementwise products of U-centered distance matrices
      dCov[i,j]=(sum(prd)-sum(diag(prd)))/(n*(n-3))    # dot-product like operator for U-centered distance matrices (page 6 of Zhu et al)
      dCov[j,i]=dCov[i,j]                              # dCov matrix is symmetric, so assign other element with same value
    }             # end loop over second dimension of X
  }               # end loop over first dimesnion of X
  
  # dCov(Y,Y)
  for (i in 1:kY)   # loop over first dimension of Y for pairs of dimensions of Y
  {
    for (j in 1:i)  # loop over second dimension of Y for pairs of dimensions of Y
    {
      prd=uc.Y[,,i]*uc.Y[,,j]                              # compute elementwise products of U-centered distance matrices
      dCov[kX+i,kX+j]=(sum(prd)-sum(diag(prd)))/(n*(n-3))  # dot-product like operator for U-centered distance matrices (page 6 of Zhu et al)
      dCov[kX+j,kX+i]=dCov[kX+i,kX+j]                      # dCov matrix is symmetric, so assign other element with same value
    }               # end loop over second dimension of X
  }                 # end loop over first dimesnion of X
  
  # dCov(X,Y)
  for (i in 1:kX)     # loop over dimensions of X
  {  
    for (j in 1:kY)   # loop over dimensions of Y
    {
      prd=uc.X[,,i]*uc.Y[,,j]                           # compute elementwise products of U-centered distance matrices
      dCov[i,kX+j]=(sum(prd)-sum(diag(prd)))/(n*(n-3))  # dot-product like operator for U-centered distance matrices (page 6 of Zhu et al)
      dCov[kX+j,i]=dCov[i,kX+j]                         # dCov matrix is symmetric, so assign other element with same value
    }                 # end loop over second dimension of X
  }                   # end loop over first dimesnion of X
  
  #################################
  # Remove large arrays that are no longer needed to save memory
  rm(uc.X)
  rm(uc.Y)
  #gc(); gc();
  
  #######################################################
  # Covert distance covariances to distance correlations
  dCor=dCov
  for (i in 1:(kX+kY))  # Loop over first matrix for pairs of all variables
  { 
    for (j in 1:i)      # Loop over second matrix for pairs of all variables
    {
      dCor[i,j]=dCov[i,j]/sqrt(dCov[i,i]*dCov[j,j])  # definition of correlation in terms of covariance and variance
      dCor[j,i]=dCor[i,j]                            # symmetry of correlation matrix
    }                   # end loop over second variable in all pairs of variables
  }                     # end loop over first variable in all pairs of variables
  
  ##################################
  # Now compute mdCov and mdCor
  const=sqrt(choose(n,2))                  # frequently used constant in Zhu et al
  odCov=sum(dCov[1:kX,kX+1:kY])*const      # defined on page 12 of Zhu et al (note 2 superscript is NOT squaring in the paper)
  odVarX=sum(dCov[1:kX,1:kX])*const        # distance-based variance of X by defintion on page 12
  odVarY=sum(dCov[kX+1:kY,kX+1:kY])*const  # distance-based variance of Y
  odCor=odCov/sqrt(odVarX*odVarY)          # definition of correlation in terms of covariance and variance
  
  #################################
  # Compute t-statistic and t-test for mdCor
  v=n*(n-3)/2                              # defined on page 18 of Zhu et al
  t.odCor=odCor/sqrt(1-odCor^2)*sqrt(v-1)  # defined on page 17 of Zhu et al
  p.odCor=2*pt(-abs(t.odCor),v-1)          # defined on page 18 and 19 of Zhu et al
  
  ################################
  # Compute t-stat and t-test for each dCor
  t.dCor=dCor/sqrt(1-dCor^2)*sqrt(v-1)    # defined on page 18 of Zhu et al
  diag(t.dCor)=NA                         # defined on page 17 of Zhu et al
  p.dCor=2*pt(-abs(t.dCor),v-1)           # defined on pages 18 and 19 of Zhu et al

  ###################################
  # result object: a list with the following components
  
  res=list(odCor=odCor,                         # overall distance correlation statistic
           t.odCor=t.odCor,                     # t-stat for overall distance correlation statistic
           p.odCor=p.odCor,                     # p-value for overall distance correlation statistic
           dCor=dCor,                           # distance-based correlation matrix for each pair of variables 
           t.dCor=t.dCor,                       # t-stat for distance-based correlation matrix
           p.dCor=p.dCor,                       # p-value for distance-based correlation matrix
           X=X,Y=Y,                             # echo input data matrices
           x.dist=x.dist,                       # echo input distance metric for X
           y.dist=y.dist)                       # echo input distance metric for Y

  class(res)="dcor"
  
  return(res)                # return results
}
