# some tools for data munging

sum1 = function(x){
  # ensure sums to 1
  return(x/sum(x))
}

cum1 = function(x){
  # cumsum spanning (0,1]
  return(cumsum(x)/sum(x))
}

int.cut = function(x,low){
  # cut with simplified labels assuming integer data
  # e.g. int.cut(seq(6),c(1,2,3,5)) -> c("1","2","3 - 4", "3 - 4","5 +","5 +")
  high = c(low[2:len(low)]-1,Inf)
  labels = gsub('-Inf','+',ifelse(low==high,low,paste0(low,'-',high)))
  x.cut = cut(x,breaks=c(low,Inf),labels=labels,right=FALSE)
}

tabulate0 = function(x,xmax=max(x)){
  c(sum(x==0),tabulate(x,xmax))
}

sample.wtd = function(x,w,n=len(x)){
  # same as sample(x,prob=w,repl=FALSE) but much faster
  # https://stackoverflow.com/a/15205104/5228288
  key = runif(len(x)) ^ (1/w)
  x.sample = x[order(key,decreasing=TRUE)][1:n]
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

list.str = function(x){
  # convert a list to strings (characters)
  # e.g. list.str(list(a=1,b=2)) -> c('a = 1','b = 2')
  sapply(names(x),function(name){
    paste(name,'=',paste(x[[name]],collapse=','))
  })
}

r.fun = function(fun,...,shift=0,rmin=NULL,rmax=NULL){
  # pre-specify some arguments to random number generating fun
  # e.g. f = r.fun(runif,min=1,max=2); f(n=10) -> runif(n=10,min=1,max=2)
  args = list(...)
  rfun = function(...){
    r = kw.call(fun,args,...) + shift
    if (!is.null(rmin)){ r = pmax(r,rmin) }
    if (!is.null(rmax)){ r = pmin(r,rmax) }
    return(r)
  }
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
