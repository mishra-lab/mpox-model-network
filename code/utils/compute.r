# some tools for (expensive) compute stuff

.n.cores = 7 # allow default global override as par.lapply is often nested within functions

par.lapply = function(...,cores=.n.cores,.par=TRUE){
  # simple wrapper for parallel::mclapply with some default arguments
  if (.par){
    parallel::mclapply(...,mc.cores=cores,mc.set.seed=FALSE)
  } else {
    lapply(...)
  }
}
