# some tools for (expensive) compute stuff

library('parallel')

.n.cores = 7 # allow default global override as par.lapply is often nested within functions

par.lapply = function(...,cores=.n.cores){
  # simple wrapper for parallel::mclapply with some default arguments
  mclapply(...,
    mc.cores    = cores,
    mc.set.seed = FALSE
  )
}

par.funs = function(fun.args,...){
  # e.g. par.funs(list(list(sum,2),list(sum,2,3)),1) -> list(3,6)
  par.lapply(fun.args,function(fa){
    do.call(fa[[1]],c(fa[1+seqn(len(fa)-1)],list(...)))
  })
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
