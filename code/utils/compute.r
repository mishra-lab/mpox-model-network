# some tools for (expensive) compute stuff

library('parallel')

.n.cores = 7 # allow default global override as par.lapply is often nested within functions

par.lapply = function(...,cores=.n.cores,.par=TRUE){
  # simple wrapper for parallel::mclapply with some default arguments
  if (.par){
    mclapply(...,mc.cores=cores,mc.set.seed=FALSE)
  } else {
    lapply(...)
  }
}

par.funs = function(fun.args,...){
  # e.g. par.funs(list(list(sum,2),list(sum,2,3)),1) -> list(3,6)
  par.lapply(fun.args,function(fa){
    do.call(fa[[1]],c(fa[1+seqn(len(fa)-1)],list(...)))
  })
}

optim.brute = function(fun,rng,n=21,d=0,.par=TRUE){
  # find x which minimizes fun within range rng, by brute force
  # evaluate fun at n points within rng
  # d > 0 recurses to improve precision by ~ 1/(n-1)
  x.i = seq(rng[1],rng[2],l=n)
  e.i = unlist(par.lapply(x.i,fun,.par=.par))
  if (d > 0){
    x.rng = sort(x.i[match(sort(e.i)[2:3],e.i)]) # use adjacent points to minimum, in case skew
    x.min = optim.brute(fun,x.rng,n=n,d=d-1,.par=.par) # recurse with new range and d-1
  } else {
    x.min = x.i[min(e.i)==e.i] # best estimate
  }
}

load.bar = function(i,N,width=100){
  if (missing(i) && missing(N)){ i=100; N=100; fin=TRUE } else {fin=FALSE}
  if (0 == (i %% (N/width))) {
    n.done = round(i/N*width)
    n.todo = width - n.done
    bksp = paste0(rep('\b',width+2),collapse='')
    done = paste0(rep('#',n.done),collapse='')
    todo = paste0(rep('-',n.todo),collapse='')
    cat(paste0(bksp,'[',done,todo,']'))
  }
  if (fin){ cat('\n') }
}

profile.code = function(expr,fname='.tmp/Rprof.out',interval=.01,...){
  Rprof(fname,interval=interval,...)
  expr;
  Rprof(NULL)
  summaryRprof(fname)
}
