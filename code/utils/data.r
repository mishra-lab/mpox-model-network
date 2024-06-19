# some tools for data munging

sum1 = function(x){
  # ensure sums to 1
  return(x/sum(x))
}

int.cut = function(x,low){
  # cut with simplified labels assuming integer data
  # e.g. int.cut(seq(6),c(1,2,3,5)) -> c("1","2","3 - 4", "3 - 4","5 +","5 +")
  high = c(low[2:len(low)]-1,Inf)
  labels = gsub('-Inf','+',ifelse(low==high,low,paste0(low,'-',high)))
  x.cut = cut(x,breaks=c(low,Inf),labels=labels,right=FALSE)
}

chunk.fun = function(x,l,n,fun){
  # apply a function over n regular chunks of length l
  # e.g. chunk.fun(1:9,3,fun=cumsum) -> c(1,3,6, 4,9,15, 7,15,24)
  if (missing(l)){ l = len(x)/n }
  if (missing(n)){ n = len(x)/l }
  g = rep(1:n,each=l)
  x.c = unlist(lapply(split(x,g),fun))
}

iter.prop.fit = function(x,xs,tol=.1,iter.max=100){
  # re-scales x row & column-wise to match the margins (row/col-sums) xs
  # also ensures the result is exactly symmetric by averaging with the transpose
  for (k in seqn(iter.max)){
    s1 = rowSums(x)
    x = x * rep(xs/s1,times=len(xs))
    s2 = colSums(x)
    x = x * rep(xs/s2,each=len(xs))
    if (all(abs(s1-xs) < tol)){ return((x+t(x))/2) }
  }
  stop('No solution to iter.prop.fit with tol = ',tol)
}

round.sum = function(x,xs=NULL){
  # rounds x, preserving xs = round(sum(x)) by rounding closest elements up / down
  # e.g. round.sum(c(3.3,3.4,3.3)) -> c(3,4,3)
  if (is.null(xs)){ xs = round(sum(x)) }
  x1 = round(x)
  xf = x %% 1
  ds = sum(x1) - xs
  if (ds < 0){
    i = head(order(xf,decreasing=TRUE),-ds)
    x1[i] = x1[i] + 1
  }
  if (ds > 0){
    xf[x1==0] = NA
    i = head(order(xf,decreasing=FALSE),ds)
    x1[i] = x1[i] - 1
  }
  return(x1)
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

sample.wtd = function(x,w,n=len(x)){
  # same as sample(x,prob=w) but much faster
  # https://stackoverflow.com/a/15205104/5228288
  key = runif(len(x)) ^ (1/w)
  x.sample = x[order(key,decreasing=TRUE)][1:n]
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

formula.fun = function(y,...){
  # construct a formula from a list of variables
  # e.g. formula fun('a','b','c') -> 'a ~ b + c'
  formula(paste(y,'~',paste(unique(c(...)),collapse=' + ')))
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
  if (is.null(m)){ m = self.name(names(x)) }
  for (name in names(x)){
    y[i %in% x[[name]]] = m[name]
  }
  return(y)
}

date.vec = function(t.vec,date.t0){
  # e.g. date.vec(1:3,'2000-01-01') -> c('2000-01-01','2000-01-02','2000-01-03')
  as.Date(date.t0) + t.vec - 1
}

interp.fun = function(x,y,pad=FALSE){
  # set-up inputs for spinterp so only xp is needed
  if (pad){
    x = c(x[1]-1,x,x[len(x)]+1)
    y = c(y[1],y,y[len(y)])
  }
  function(xp){ pracma::spinterp(x,y,xp) }
}

gsmooth = function(x,bw=1){
  # gaussian smoothing
  m = round(4*bw)
  kern = sum1(dnorm(-m:+m,0,bw))
  xs = convolve(x,kern,type='open')[(m+1):(m+len(x))]
}

t.doubling = function(x,lag=1,bw=1,clip=FALSE,pad=TRUE){
  # estimate doubling time of x after gaussian smoothing
  t2x = lag/diff(log(gsmooth(x,bw),base=2),lag=lag)
  if (clip){ t2x = pmin(clip,pmax(-clip,t2x)) }
  if (pad){ t2x = c(t2x,rep(t2x[len(t2x)],lag)) }
}
