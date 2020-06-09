surv.dist <-
function(stime.evnt) # data.frame with survival time in column 1 and event status in column 2(0 indicates censored, 1 indicates event)
{
  n=nrow(stime.evnt)                 # sample size
  res=matrix(0,n,n)                  # initialize matrix
  evnt.indx=which(stime.evnt[,2]>0)  # index of rows with an event
  
  time.pts=sort(unique(stime.evnt[evnt.indx,1])) # unique event times
  for (i in 1:length(time.pts)) # loop over unique event times
  {
    early.event=(stime.evnt[,1]<=time.pts[i])&
      (stime.evnt[,2]==1)                       # subjects with event before this event time
    
    late.event1=(stime.evnt[,1]>time.pts[i])    # subjects with longer follow-up than this event time
    
    late.event2=(stime.evnt[,1]==time.pts[i])&  # subjects censored at this event time
      (stime.evnt[,2]==0)
    
    late.event=(late.event1|late.event2)        # subjects that definitely have event after this event time
    
    temp=matrix(0,n,n)                          # initialize distance coefficient matrix to all zeroes
    temp[early.event,early.event]=0             # assign distance coefficient matrix to -1 for pairs of subjects with earlier events
    temp[late.event,late.event]=0               # assign distance coefficient matrix to -1 for pairs of subjects with later events
    temp[early.event,late.event]=1              # assign distance coefficient matrix to +1 for pairs of subjects with early and late event
    temp[late.event,early.event]=1              # assign distance coefficient matrix to +1 for pairs of subjects with early and late event
    
    res=res+temp                # add result for this time point to cumulative coefficient matrix
  }                             # end loop over risk sets
  
  res=as.dist(res)
  return(res) # return the matrix
}
