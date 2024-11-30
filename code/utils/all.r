# -----------------------------------------------------------------------------
# options

options(
  stringsAsFactors=FALSE,
  showNCalls=500,
  nwarnings=1e4,
  scipen=99,
  width=200)

# -----------------------------------------------------------------------------
# files

root.path = function(...,create=FALSE){
  # make a file path starting one level above /code/
  # e.g. root.path('foo','bar') -> '.../{projectroot}/foo/bar'
  root = strsplit(file.path(getwd(),''),file.path('','code',''))[[1]][1]
  path = file.path(root,...)
  if (create & !dir.exists(dirname(path))){ dir.create(dirname(path),recursive=TRUE) }
  return(path)
}

tmp.out = function(out,slug,path='.tmp',ext='.out'){
  # print some text (out) and also save it to a temporary file; also, path can be a vector
  # e.g. path=c('a','b') -> 'a/b/{slug}-1970-01-01.out'
  file.out = kw.call(file.path,c(as.list(path),paste0(slug,'-',Sys.Date(),ext)))
  print(out); sink(file.out); print(out); sink()
}

# -----------------------------------------------------------------------------
# plotting

suppressPackageStartupMessages({
  library('ggplot2')
  library('reshape2')
})

.plot.ext = 'pdf'

fig.save = function(...,g=last_plot(),w=7,h=7,ext=.plot.ext){
  # wrapper for ggsave
  fname = paste0(...,'.',ext)
  print(paste('saving:',fname))
  ggsave(fname,plot=g,w=w,h=h)
}

# -----------------------------------------------------------------------------
# misc

len = length
lens = lengths
seqn = seq_len

sum1 = function(x){ x/sum(x) }

cum1 = function(x){ cumsum(x)/sum(x) }

tabulate0 = function(x,...){ c(sum(x==0),tabulate(x,...)) }

dn.array = function(dn,x=NA){
  # easily construct array from dimnames
  A = array(x,dim=sapply(dn,len),dimnames=dn)
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

row.select = function(x,select=list(),...,.return='x'){
  # e.g. row.select(data.frame(x=c(1,2,3)),list(x=c(1,3))) -> data.frame(x=c(1,3))
  b = rep(TRUE,nrow(x))
  select = list.update(select,...)
  for (name in names(select)){
    b = b & x[[name]] %in% select[[name]]
  }
  return(switch(.return,
    'x' = x[b,],
    'b' = b))
}

split.col = function(x,col,new.names,del=TRUE){
  # e.g. split.col(data.frame(y=c('a.1','a.2','b.1','b.2')),'y',c('char','num'))
  #   -> data.frame(char=c('a','a','b','b'),num=c('1','2','1','2'))
  new.cols = as.data.frame(do.call(rbind,strsplit(as.character(x[[col]]),'\\.')))
  names(new.cols) = new.names
  x = cbind(x,new.cols)
  if (del){ x[[col]] = NULL }
  return(x)
}

lookup.map = function(i,x,m=NULL){
  # e.g. lookup.map(1:5,list(a=1:2,b=3:5)) -> c('a','a','b','b','b')
  # e.g. lookup.map(1:5,list(a=1:2,b=3:5),list(a='A',b='B')) -> c('A','A','B','B','B')
  y = NA * i
  if (is.null(m)){ m = setNames(names(x),names(x)) }
  for (name in names(x)){
    y[i %in% x[[name]]] = m[name]
  }
  return(y)
}

sample.wtd = function(x,w,n=len(x)){
  # same as sample(x,prob=w,repl=FALSE) but much faster
  # https://stackoverflow.com/a/15205104/5228288
  key = runif(len(x)) ^ (1/w)
  x.sample = x[order(key,decreasing=TRUE)][1:n]
}

r.fun = function(fun,args=list(),...,shift=0,rmin=NULL,rmax=NULL){
  # pre-specify some arguments to random number generating fun
  # e.g. f = r.fun(runif,min=1,max=2); f(n=10) -> runif(n=10,min=1,max=2)
  args = list.update(args,...)
  rfun = function(...){
    r = kw.call(fun,args,...) + shift
    if (!is.null(rmin)){ r = pmax(r,rmin) }
    if (!is.null(rmax)){ r = pmin(r,rmax) }
    return(r)
  }
}

int.cut = function(x,low){
  # cut with simplified labels assuming integer data
  # e.g. int.cut(seq(6),c(1,2,3,5)) -> c("1","2","3 - 4", "3 - 4","5 +","5 +")
  high = c(low[2:len(low)]-1,Inf)
  labels = gsub('-Inf','+',ifelse(low==high,low,paste0(low,'-',high)))
  x.cut = cut(x,breaks=c(low,Inf),labels=labels,right=FALSE)
}

# -----------------------------------------------------------------------------
# function calls

kw.call = function(fun,kwds,...){
  # e.g. kw.call(sum,list(1,2),3) -> 6
  do.call(fun,c(list(...),kwds))
}

rbind.lapply = function(...,.par=FALSE){
  # convenience wrapper as shown below
  do.call(rbind,par.lapply(...,.par=.par))
}

.n.cores = 7 # allow default global override as par.lapply is often nested within functions

par.lapply = function(...,cores=.n.cores,.par=TRUE){
  # simple wrapper for parallel::mclapply with some default arguments
  if (.par){
    parallel::mclapply(...,mc.cores=cores,mc.set.seed=FALSE)
  } else {
    lapply(...)
  }
}
