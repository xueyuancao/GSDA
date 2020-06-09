U.center <-
function(d)
  
{
  if (nrow(d)!=ncol(d))                     # check that the matrix is square
    stop("d must be a square matrix.")
  if (any(d<0))                             # check that all elements are non-negative
    stop("all elements of d must be >=0.")
  n=nrow(d)                                 # dimension of the matrix
  rsm=rowSums(d)                            # row sums of distance matrix
  csm=colSums(d)                            # column sum of distance matrix
  osm=sum(d)-sum(diag(d))                   # overall sum of distance matrix
  rsm=matrix(rsm,n,n,byrow=F)               # row sums as a matrix
  csm=matrix(csm,n,n,byrow=T)               # column sums as a matrix
  d=d-rsm/(n-2)-csm/(n-2)+osm/((n-1)*(n-2)) # U-centering formula on page 6 of arXiv 1902.03291 paper
  return(d)
}
