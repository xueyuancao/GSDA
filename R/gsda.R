gsda <-
function(omic.data,  # matrix of omic data with variables as rows and subjects as columns
              clin.data,  # data.frame of clinical data matrix with subjects as rows and variables as columns, must include a column named "ID" with identifiers matching column names of omic.data
              vset.data,  # data.frame or matrix of gene-set assignments, one row per assignment of variable to variable set, columns "vID" for variable ID and "vset" for variable set
              clin.vars,  # vector name or numeric index of endpoint/treatment variable(s) in clin.data
              omic.dist,  # distance metric for omic data, may be "oe" (overall Euclidean), "me" (marginal Euclidean), "om" (overall Manhattan), or "mm" (marginal Manhattan)
              clin.dist)  # distance metric for clinical data, may be "oe", "me", "om", or "mm" as defined above, or "ct" (categorical) or "st" (survival time)
  
{
  # prepare data for analysis
  t0=proc.time()
  gsda.data=prep.gsda(omic.data,
                      clin.data,
                      vset.data)
  
  omic.data=gsda.data$omic.data
  clin.data=gsda.data$clin.data
  vset.data=gsda.data$vset.data
  vset.index=gsda.data$vset.index
  t1=proc.time()
  prep.time=(t1-t0)[3]
  comp.time=(t1-t1)[3]
  
  n=nrow(clin.data)          # number of subjects
  nset=nrow(vset.index)      # number of gene-sets
  p.vset=rep(0,nset)         # initialize p-value vector
  vset.list=rep("NA",nset)   # initialize vector of variable set lists 
  dCor.stats=rep(0,nset)    # initialize vector of marginal distance correlation statistics
  comp.time=rep(NA,nset)     # initialize vecto of compute time for each variable set
  
  # Loop over gene-sets
  for (i in 1:nset)
  {
    # Get variables in this variable set
    t0=proc.time()
    vset.ind=(vset.index$start.index[i]:vset.index$end.index[i]) # rows with variables in this variable set
    vset.vIDs=vset.data$vID[vset.ind]                            # ids of variables in this variable set
    vset.list[i]=paste(sort(vset.vIDs),collapse=" | ")           # string listing all variables in this variable set
    
    # Extract data matrix for this variable set
    vset.mtx=matrix(omic.data[vset.vIDs,],                       # data matrix for this variable set
                    length(vset.ind),n)    
    
    # Compute distances for this variable set
    vset.mtx=t(vset.mtx)                                         # transpose the matrix for use in dist.corr
    t1=proc.time()                                               # check clock to calculate computing times
    prep.time=prep.time+(t1-t0)[3]                               # calculate data prep time for this variable set
    
    temp.res=dist.corr(vset.mtx,                                 # compute the distance correlation stats for this variable set
                       clin.data[,clin.vars],
                       omic.dist,clin.dist)
    
    p.vset[i]=temp.res$p.odCor                                   # get the overall distance correlation p-value
    dCor.stats[i]=temp.res$odCor                                # get the overall distance correlation statistic
    t2=proc.time()                                               # get time for calculating compute time
    comp.time[i]=(t2-t1)[3]     
  } # end loop over gene-sets
  
  ################
  # result object: a data.frame with one row per gene-set and the following colums:
  res=cbind.data.frame(vset=vset.index$value,  # name of variable set (gene-set)
                       vIDs=vset.list,         # character string with list of variables in the variable set (gene-set)
                       dCor=dCor.stats,      # distance statistic for the variable set
                       p.vset=p.vset,          # p-value
                       comp.time=comp.time)    # computing time

  class(res)="gsda.result" # did not work as intended in the print function below
  
  
  return(res)
}

order.index.dset <-
function(dset,  # data.frame
                          vr)    # name or numeric index of variable by which to order and index dset
  
{
  v=dset[,vr]     # extract the variable
  ord=order(v)    # vector of indices to order by that variable
  dset=dset[ord,] # order dset by that variable
  m=nrow(dset)    # get size of data set
  new.val=which(dset[-1,vr]!=dset[-m,vr]) # find which rows have a new value of variable
  ind1=c(1,new.val+1)                     # first row with a unique value of the variable
  ind2=c(new.val,m)                       # last row with a unique value of the variable
  
  # data.frame with indices of dset that have unique values of the variable vr
  dset.ind=cbind.data.frame(start.index=ind1, 
                            end.index=ind2,
                            value=dset[ind1,vr])
  # create a list to return
  res=list(dset=dset,
           indx=dset.ind)
  # return that list
  return(res)
}

match.order.ids <-
function(id1,id2)
  
{
  temp1=cbind.data.frame(id=id1,ind1=1:length(id1)) # create data.frame with indices for id1
  temp2=cbind.data.frame(id=id2,ind2=1:length(id2)) # create data.frame with indices for id2
  res=merge(temp1,temp2,by=1)                       # merged data.frame
  
  # Warn user if some elements are not present in combined data.frame
  if (nrow(res)<length(id1))
    warning(paste0("Failed to merge ",length(id1)-nrow(res)," elements of id1."))
  
  if (nrow(res)<length(id2))
    warning(paste0("Failed to merge ",length(id2)-nrow(res)," elements of id2."))
  
  # return the result
  return(res)
}

dist.metric.name <-
function(dm)
{
  if (tolower(dm)=="oe") return("Overall Euclidean")
  if (tolower(dm)=="me") return("Marginal Euclidean")
  if (tolower(dm)=="om") return("Overall Manhattan")
  if (tolower(dm)=="mm") return("Marginal Manhattan")
  if (tolower(dm)=="ct") return("Categorical Distance")
  if (tolower(dm)=="st") return("Survival Time Distance")
}