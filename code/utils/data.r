# some tools for data munging

seq.group = function(N.i,N.g){
  # balanced rep(1:N.g) to total length N.i
  # e.g. seq.group(10,3) -> c(1,1,1,1, 2,2,2, 3,3,3)
  m = N.i/N.g
  r.0 = N.i %% N.g
  r.1 = N.g - r.0
  rep(1:N.g,times=c(rep(ceiling(m),r.0),rep(floor(m),r.1)))
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

q3.aggr = function(y,by,data,intervals=.9){
  # aggregate data[[y]] by data[[by]] to obtain median + some quantile edges
  # e.g. y='y' & intervals = c(.5,.95) yields columns y.0.5, y.0.025, y.0.25, y.0.75, y.0.975
  q.vec = c(.5, rev(1-intervals)/2, 1-(1-intervals)/2)
  f = formula(paste(y,'~',paste(by,collapse='+')))
  data.aggr = aggregate(f,data,quantile,q.vec)
  y.aggr = as.data.frame(data.aggr[[y]])
  data.aggr[[y]] = NULL
  names(y.aggr) = paste0(y,'.',q.vec)
  return(cbind(data.aggr,y.aggr))
}

q.interval = function(interval){
  # e.g. interval = .95 -> c(0.025, 0.975)
  ci = c((1-interval)/2, 1-(1-interval)/2)
}

