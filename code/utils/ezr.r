# some tools to make life easier

len = length

seqn = seq_len

sum1 = function(x){
  # ensure sums to 1
  return(x/sum(x))
}

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

dn.array = function(dn,x=NA){
  # easily construct array from dimnames
  A = array(x,dim=sapply(dn,len),dimnames=dn)
}
