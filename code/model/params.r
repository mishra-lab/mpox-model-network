
def.params = function(seed=NULL,...){
  set.seed(seed)
  P = list()
  P$seed           = seed
  P$lab.city       = c('A','B','C') # city labels
  P$N.city         = c(6,5,4) * 50 # city pop sizes
  P$N              = sum(P$N.city) # total pop size
  P$net.type       = 'fit'
  P$net.args.city  = list(
    'A' = def.params.msm.net(P$N.city[1]),
    'B' = def.params.msm.net(P$N.city[2]),
    'C' = def.params.msm.net(P$N.city[3]))
  P$e.bridge       = .01 # fraction of bridges (contacts) between cities
  P$mix.bridge     = matrix(c(0,3,1, 3,0,1, 1,1,0),c(3,3)) # relative distribution of bridges
  P$N.I0.city      = c(10,0,0) # number initially infected per-city
  P$dur.exp        = 7 # average duration in exposed
  P$dur.inf        = 7 # average duration infectious (considering isolation)
  P$beta           = .9 # overall probability of transmission
  P$net.dur        = 30 # TODO: should be 180 per data, but haven't modelled repeat contacts ...
  P$vax.eff.dose   = c(.60,.90) # vaccine effectiveness by dose
  N.vax = .10 * P$N
  P$vax.args.phase = list( # vaccination config
    '1'  = list(t0=30,dur=20,N=N.vax,dose=1,w.city=sum1(P$N.city),w.attr=NULL),
    '2a' = list(t0=60,dur=20,N=N.vax,dose=2,w.city=sum1(P$N.city),w.attr=NULL),
    '2b' = list(t0=60,dur=20,N=N.vax,dose=1,w.city=sum1(P$N.city),w.attr=NULL))
  P = list.update(P,...) # override any of the above
  if (is.null(P$G)){
    P$G = gen.multi.city(P) # generate network (expensive)
  }
  return(P)
}

def.params.s = function(seeds){
  # run def.params for a number (or vector) of seeds, parallel because net gen is expensive
  if (len(seeds)==1){ seeds = seqn(seeds) }
  return(par.lapply(seeds,def.params))
}

def.params.msm.net = function(N){
  g.par = list(shape=0.255,rate=0.032,shift=0.913) # empiric
  P = list()
  P$deg.mean = 3.412 + 0.323 * max(0,log10(N)) # empiric
  P$fitness  = rgamma(N,shape=g.par$shape,rate=g.par$rate) + g.par$shift
  P$adj.power = -.628 # empiric
  return(P)
}
