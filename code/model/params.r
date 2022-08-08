
def.params = function(seed=NULL,N=1000,...){
  set.seed(seed)
  P = list()
  # independent parameters (mostly)
  P$seed           = seed
  P$N              = N # pop size total
  P$net.dur        = 6*30 # period of time reflected in the sexual network
  P$net.params     = def.params.net(P$N) # params for network
  P$N.I0           = 10 # number initially infected
  P$dur.exp        = 8.5 # average duration in exposed (latent period) - expon distrib
  P$dur.inf        = 7 # average duration infectious (considering isolation) - expon distrib
  P$beta           = .90 # probability of transmission (per contact)
  P$vax.eff.dose   = c(.60,.90) # vaccine effectiveness by dose
  P$N.V0           = c(.00,.00) * P$N # total number initially vaccinated by dose
  P = list.update(P,...) # override any of the above
  # conditional parameters
  if (is.null(P$G)){ P$G = make.net(P$net.params) } # generate the sexual network
  P$beta.health = c('S'=1,'E'=0,'I'=0,'R'=0,'V1'=1-P$vax.eff.dose[1],'V2'=1-P$vax.eff.dose[2])
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
  P.net$par.gam.shape = 0.255 # gamma distrib: partners in 6 months
  P.net$par.gam.rate  = 0.032 # gamma distrib: partners in 6 months
  P.net$par.gam.shift = 0.913 # gamma distrib: partners in 6 months
  P.net$main.i.frac = .2 # fraction of pop who have main partners
  P.net$main.w.par.power = -1 # when chosing who has main partners: weights = p6m ^ power
  P.net$main.sex.gam.shape = 5  # gamma distrib: sex per main partner in 6 months
  P.net$main.sex.gam.rate  = .2 # gamma distrib: sex per main partner in 6 months
  P.net$main.sex.gam.shift = 1  # gamma distrib: sex per main partner in 6 months
  P.net$casu.sex.gam.shape = .3 # gamma distrib: pex per casu partner in 6 months
  P.net$casu.sex.gam.rate  = .3 # gamma distrib: pex per casu partner in 6 months
  P.net$casu.sex.gam.shift = 1  # gamma distrib: pex per casu partner in 6 months
  return(P.net)
}

make.net = function(P.net){
  # the sexual network reflects all contacts (sex) occuring in P$net.dur days (6 months)
  # including multiple contacts per partnership
  i = seqn(P.net$N)
  # sample total partners in 6 months
  par.i = round(rgamma(P.net$N,P.net$par.gam.shape,P.net$par.gam.rate) + P.net$par.gam.shift)
  par.i = degrees.balanced(par.i)
  # generate main partners
  i.main = sample(i,round(P.net$main.i.frac*P.net$N),p=par.i^P.net$main.w.par.power) # who
  ii.e.main = edges.random(i.main,shuffle=FALSE) # edges = pairs
  sex.e.main = round(rgamma(nrow(ii.e.main),P.net$main.sex.gam.shape,P.net$main.sex.gam.rate) + P.net$main.sex.gam.shift)
  ii.e.main.sex = edges.repeated(ii.e.main,sex.e.main) # duplicate edges (sex) within each pair
  # generate casual partnerss
  par.i.casu = par.i
  par.i.casu[i.main] = par.i.casu[i.main] - 1 # non-main partners
  ii.e.casu = edges.unloop(edges.from.degrees(i,par.i.casu)) # edges = pairs
  sex.e.casu = round(rgamma(nrow(ii.e.casu),P.net$casu.sex.gam.shape,P.net$casu.sex.gam.rate) + P.net$casu.sex.gam.shift)
  ii.e.casu.sex = edges.repeated(ii.e.casu,sex.e.casu) # duplicate edges (sex) within each pair
  # all contacts
  ii.e = rbind(ii.e.main.sex,ii.e.casu.sex)
  # attributes
  g.attr = list()
  i.attr = list()
  i.attr$par = par.i
  e.attr = list()
  if (.debug){ # expensive / not required
    i.attr$sex = degrees.from.edges(i,ii.e)
    i.attr$main = factor(i %in% i.main,c(T,F),c('Yes','No'))
    e.attr$type = factor(c(rep('main',nrow(ii.e.main.sex)),rep('casu',nrow(ii.e.casu.sex))))
    kr.e = index.repeated.edges(ii.e)
    e.attr$sex.index = kr.e[,1]
    e.attr$sex.total = kr.e[,2]
  }
  # graph object
  # note: G$deg.i is misused here: should be total sex, but we store total partners
  G = graph.obj(ii.e=ii.e,i=i,deg.i=par.i,g.attr=g.attr,i.attr=i.attr,e.attr=e.attr)
  if (.debug){ G$attr$g$layout = graph.layout(G) } # pre-compute consistent layout if needed
  return(G)
}