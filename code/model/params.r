
def.params = function(seed=NULL,...){
  set.seed(seed)
  P = list()
  P$seed           = seed
  P$N.city         = c(10000,0,0) # pop sizes by city
  P$N              = sum(P$N.city) # pop size total
  P$net.dur        = 6*30 # period of time reflected in the sexual network
  P$net.params.city  = list( # params for networks in each city
    'A' = def.params.net.city(P$N.city[1]),
    'B' = def.params.net.city(P$N.city[2]),
    'C' = def.params.net.city(P$N.city[3]))
  # TODO: rebuild bridges, later
  # P$e.bridge       = .00 # fraction of bridges (contacts) between cities
  # P$mix.bridge     = matrix(c(0,3,1, 3,0,1, 1,1,0),c(3,3)) # relative distribution of bridges
  P$N.I0.city      = c(10,0,0) # number initially infected per-city
  P$dur.exp        = 8.5 # average duration in exposed (latent period) - expon distrib
  P$dur.inf        = 7 # average duration infectious (considering isolation) - expon distrib
  P$beta           = .90 # probability of transmission (per contact)
  P$vax.eff.dose   = c(.60,.90) # vaccine effectiveness by dose
  N.vax = .00 * P$N # total number vaccinated
  P$vax.params.phase = list( # vaccination config
    '01' = def.params.vax.phase(dose=1,t=0,N.total=1*N.vax,w.city=sum1(P$N.city),w.attr=NULL),
    '02' = def.params.vax.phase(dose=2,t=0,N.total=0*N.vax,w.city=sum1(P$N.city),w.attr=NULL))
  P = list.update(P,...) # override any of the above
  if (is.null(P$G)){ # generate the sexual network
    P$G = make.net.city(P$net.params.city$A,'A') # TEMP
    # P$G = make.net.multi.city(P) # TODO: rebuild
  }
  P$seed.state = .Random.seed # current state
  return(P)
}

def.params.s = function(seeds){
  # run def.params for a number (or vector) of seeds, parallel because net gen is expensive
  if (len(seeds)==1){ seeds = seqn(seeds) }
  return(par.lapply(seeds,def.params))
}

def.params.vax.phase = function(dose,t,N.total,w.city,w.attr){
  P.vax = list()
  P.vax$dose = dose # 1 or 2
  P.vax$t = t # days of the vaccination campaign, e.g. 31:40
  P.vax$N.total = N.total # total vaccines in this phase
  P.vax$N.day.city = do.call(cbind,lapply(w.city*N.total,groups.even,N.g=len(t))) # doses / day, city
  P.vax$w.city = w.city # allocation by city (exact)
  P.vax$w.attr = w.attr # allocation by individual attributes (random)
  return(P.vax)
}

def.params.net.city = function(N){
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

make.net.city = function(P.net,city){
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
  g.attr$city = city
  i.attr = list()
  i.attr$city = rep(city,P.net$N)
  i.attr$par = par.i
  e.attr = list()
  if (.debug){ # expensive / not required
    i.attr$sex = degrees.from.edges(i,ii.e)
    i.attr$main = factor(i %in% i.main,c(T,F),c('Yes','No'))
    e.attr$type = factor(c(rep('main',nrow(ii.e.main.sex)),rep('casu',nrow(ii.e.casu.sex))))
    e.attr$city = rep(city,nrow(ii.e))
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