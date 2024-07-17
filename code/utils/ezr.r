# some tools to make life easier

len = length
lens = lengths
seqn = seq_len

kw.call = function(fun,kwds,...){
  # e.g. kw.call(sum,list(1,2),3) -> 6
  do.call(fun,c(list(...),kwds))
}

list.update = function(x,xu=list(),...){
  # e.g. list.update(list(a=1,b=2),xu=list(a=3),b=4) -> list(a=3,b=4)
  args = c(xu,list(...))
  for (name in names(args)){
    x[[name]] = args[[name]]
  }
  return(x)
}

rename.cols = function(x,...){
  # e.g. rename.cols(X,a='apples') would rename column 'a' to 'apples'
  args = list(...)
  for (name in names(args)){
    x[[args[[name]]]] = x[[name]]
    x[[name]] = NULL
  }
  return(x)
}

remove.cols = function(x,...){
  # e.g. remove.cols(X,'a','b') would remove columns 'a' and 'b'
  args = c(...)
  for (name in args){ x[[name]] = NULL }
  return(x)
}

rbind.lapply = function(...,.par=FALSE){
  # convenience wrapper as shown below
  do.call(rbind,par.lapply(...,.par=.par))
}

dn.array = function(dn,x=NA){
  # easily construct array from dimnames
  A = array(x,dim=sapply(dn,len),dimnames=dn)
}
