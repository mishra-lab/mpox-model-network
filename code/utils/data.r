# some tools for data munging

int.cut = function(x,low){
  # cut with simplified labels assuming integer data
  # e.g. int.cut(seq(6),c(1,2,3,5)) -> c("1","2","3 - 4", "3 - 4","5 +","5 +")
  high = c(low[2:len(low)]-1,Inf)
  labels = gsub('- Inf','+',ifelse(low==high,low,paste(low,'-',high)))
  x.cut = cut(x,breaks=c(low,Inf),labels=labels,right=FALSE)
}

cumsum.len = function(x,l,n){
  if (missing(l)){ l = len(x)/n }
  g = rep(1:n,each=l)
  x.cs = unlist(lapply(split(x,g),cumsum))
}

groups.even = function(N.i,N.g){
  # balanced rep(1:N.g) to total length N.i
  # e.g. seq.group(10,3) -> c(1,1,1,1, 2,2,2, 3,3,3)
  m = N.i/N.g
  r.0 = N.i %% N.g
  r.1 = N.g - r.0
  N.i.g = c(rep(ceiling(m),r.0),rep(floor(m),r.1))
}

sample.i = function(i,p){
  i[runif(len(i)) < p]
}

sample.strat = function(x,n,strat=NULL,weights=NULL,cap.n=TRUE){
  # sampling, possibly with stratification & weights; also maybe ensure n.i < len(x.i)
  if (is.null(weights)){
    weights = rep(1,len(x))
  }
  if (is.null(strat)){
    x.sample = sample(x,n,p=weights)
  } else {
    x.sample = unname(unlist(mapply(function(x.i,n.i,w.i){
      if (cap.n){ n.i = min(n.i,len(x.i)) }
      x.i.sample = sample(x=x.i,size=n.i,p=w.i)
    },split(x,strat),n,split(weights,strat))))
  }
}

row.select = function(x,select=list(),...){
  # e.g. row.select(data.frame(x=c(1,2,3)),list(x=c(1,3))) -> data.frame(x=c(1,3))
  i = rep(TRUE,nrow(x))
  select = list.update(select,...)
  for (name in names(select)){
    i = i & x[[name]] %in% select[[name]]
  }
  return(x[i,])
}

join.str = function(...){
  # e.g. join.str('a','b') -> 'a.b' -- to match split.col
  paste(...,sep='.')
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
