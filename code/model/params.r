
def.params = function(seed=NULL,N=1000,...){
  set.seed(seed)
  P = list()
  # independent parameters (mostly)
  P$seed           = seed
  P$N              = N # pop size total
  P$net.params     = def.params.net(P$N) # params for network
  P$N.I0           = 10 # number initially infected
  P$dur.EI.rfun   = r.fun(rlnorm,meanlog=2.09,sdlog=0.46,rmin=3,rmax=21) # incubation period
  P$dur.IR.rfun   = r.fun(rgamma,shape=36,scale=0.58,rmin=14,rmax=28) # infectious period
  P$dur.IH.rfun   = r.fun(rgamma,shape=1.23,scale=4.05,rmin=2,rmax=20) # non-isolated period
  P$p.asymp        = .15 # proportion of cases asymptomatic (never isolate)
  P$beta           = .90 # probability of transmission (per contact)
  P$vax.eff.dose   = c(.85,.88) # vaccine effectiveness by dose
  P$N.V0           = c(.00,.00) * P$N # total number initially vaccinated by dose
  P$p.detect.t     = interp.fun(c(0,30),c(0,.85),pad=TRUE) # probability of detection vs t
  P = list.update(P,...) # override any of the above
  # conditional parameters
  if (is.null(P$G)){ P$G = make.net(P$net.params) } # generate the sexual network
  P$beta.health = P$beta * # transmission prob by health state (susceptibility)
    c('S'=1,'E'=0,'I'=0,'H'=0,'R'=0,'V1'=1-P$vax.eff.dose[1],'V2'=1-P$vax.eff.dose[2])
  P$seed.state = .Random.seed # current state
  return(P)
}

def.params.s = function(seeds,...){
  # run def.params for a number (or vector) of seeds, parallel because net gen is expensive
  if (len(seeds)==1){ seeds = seqn(seeds) }
  P.s = par.lapply(seeds,def.params,...)
}

def.params.net = function(N){
  P.net = list()
  P.net$N = N
  P.net$dur = 6*30 # period of time reflected in the sexual network
  # will revise below based on p6m partners, stratified by had vs did not have main
  P.net$par.rfun = r.fun(rgamma,shape=0.255,rate=0.032,shift=0.913) # partners in p6m
  P.net$main.i.frac = .2 # fraction of pop who have main partners
  P.net$main.w.par.power = -1 # when chosing who has main partners: weights = p6m ^ power
  P.net$main.sex.rfun = r.fun(rgamma,shape=5,rate=.2,shift=1) # sex per main partner in p6m
  P.net$casu.sex.rfun = r.fun(rgamma,shape=.3,rate=.3,shift=1) # sex per main partner in p6m
  return(P.net)
}

make.net = function(P.net){
  # the sexual network reflects all contacts (sex) occuring in P$net.dur days (6 months)
  # including multiple contacts per partnership
  i = seqn(P.net$N)
  # sample total partners in 6 months
  par.i = round(P.net$par.rfun(P.net$N))
  par.i = degrees.balanced(par.i)
  # generate main partners
  i.main = sample(i,round(P.net$main.i.frac*P.net$N),p=par.i^P.net$main.w.par.power) # who has main
  ii.e.main = edges.random(i.main,shuffle=FALSE) # edges = pairs
  sex.e.main = round(P.net$main.sex.rfun(nrow(ii.e.main)))
  # generate casual partners
  par.i.casu = par.i
  par.i.casu[i.main] = par.i.casu[i.main] - 1 # non-main partners
  ii.e.casu = edges.unloop(edges.from.degrees(i,par.i.casu)) # edges = pairs
  sex.e.casu = round(P.net$casu.sex.rfun(nrow(ii.e.casu)))
  # all contacts
  ii.e = rbind(ii.e.main,ii.e.casu)
  sex.e = c(sex.e.main,sex.e.casu)
  # attributes
  g.attr = list()
  g.attr$dur = P.net$dur
  i.attr = list()
  i.attr$par = par.i
  e.attr = list()
  e.attr$sex = sex.e
  if (.debug){ # expensive / not required
    # i.attr$sex = degrees.from.edges(i,ii.e) # TODO: re-implement using sex.e
    i.attr$main = factor(i %in% i.main,c(T,F),c('Yes','No'))
    e.attr$type = factor(c(rep('main',nrow(ii.e.main)),rep('casu',nrow(ii.e.casu))))
  }
  # graph object
  G = graph.obj(ii.e=ii.e,i=i,deg.i=par.i,g.attr=g.attr,i.attr=i.attr,e.attr=e.attr)
  # TODO: this results in .Random.seed depends on .debug: maybe move this after .Random.seed saved
  if (.debug){ G$attr$g$layout = graph.layout.fr(G) } # pre-compute consistent layout if needed
  return(G)
}