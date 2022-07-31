

def.params = function(seed=NULL,...){
  set.seed(seed)
  P = list()
  P$seed           = seed
  P$lab.city       = c('A','B','C') # city labels
  P$N.city         = c(6,5,4) * 50 # city pop sizes
  P$N              = sum(P$N.city) # total pop size
  P$net.type       = 'ff' # network generation type for each city, and per-city args (below)
  P$net.args.city  = list('A'=list(fw=.18,bw=2),'C'=list(fw=.20,bw=2),'B'=list(fw=.22,bw=2))
  P$e.bridge       = .01 # fraction of bridges (contacts) between cities
  P$mix.bridge     = matrix(c(0,3,1, 3,0,1, 1,1,0),c(3,3)) # relative distribution of bridges
  P$N.I0.city      = c(10,0,0) # number initially infected per-city
  # P$I0.city        = c(10,0,0) # per-city
  P$Z0.max         = 7 # Z for I0 sampled from runif(1,Z0.max)
  P$inf.dur        = 28 # duration of infectiousness
  P$inf.pz         = sum1(dgamma(seq(P$inf.dur),13,1)) # p-trans per-day post exposure
  P$beta           = .5 # overall probability of transmission
  P$vax.t0         = 30 # start vaccination roll-out
  P$vax.dur        = 10 # duration of vaccination roll-out
  P$N.vax.city.day = round(P$N.city * c(.40,.20,.00) / P$vax.dur) # number vax per-city per day
  P$w.vax.attr     = 'degree' # attribute by which to weight vaccine allocation
  P = list.update(P,...) # override any of the above
  if (is.null(P$G)){
    P$G = gen.multi.city(P) # generate network (expensive)
  }
  return(P)
}

def.params.s = function(seeds){
  # run def.params for a number (or vector) of seeds, parallel because net gen is expensive
  if (len(seeds)==1){ seeds = seq(seeds) }
  return(par.lapply(seeds,def.params))
}