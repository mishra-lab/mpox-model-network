
def.params = function(seed=NULL,...){
  set.seed(seed)
  P = list()
  P$seed           = seed
  P$N.city         = c(1000,0,0) # city pop sizes
  P$N              = sum(P$N.city) # total pop size
  P$net.dur        = 180 # period of time reflected in the net
  P$net.params.city  = list(
    'A' = def.params.net.city(P$N.city[1]),
    'B' = def.params.net.city(P$N.city[2]),
    'C' = def.params.net.city(P$N.city[3]))
  # TODO: rebuild bridges, later
  # P$e.bridge       = .00 # fraction of bridges (contacts) between cities
  # P$mix.bridge     = matrix(c(0,3,1, 3,0,1, 1,1,0),c(3,3)) # relative distribution of bridges
  P$N.I0.city      = c(10,0,0) # number initially infected per-city
  P$dur.exp        = 7 # average duration in exposed
  P$dur.inf        = 7 # average duration infectious (considering isolation)
  P$beta           = .9 # overall probability of transmission
  P$vax.eff.dose   = c(.60,.90) # vaccine effectiveness by dose
  N.vax = .10 * P$N
  P$vax.params.phase = list( # vaccination config
    '01' = def.params.vax.phase(dose=1,t=0,N.total=1*N.vax,w.city=sum1(P$N.city),w.attr=NULL),
    '02' = def.params.vax.phase(dose=2,t=0,N.total=0*N.vax,w.city=sum1(P$N.city),w.attr=NULL))
  P = list.update(P,...) # override any of the above
  if (is.null(P$G)){
    P$G = make.net.city(P$net.params.city$A,'A') # TEMP
    # P$G = make.net.multi.city(P) # TODO: rebuild
  }
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
  P.vax$t = t # days of the vaccination campaign
  P.vax$N.total = N.total # total vaccines in this phase
  P.vax$N.day.city = do.call(cbind,lapply(w.city*N.total,groups.even,N.g=len(t)))
  P.vax$w.city = w.city # allocation by city (exact)
  P.vax$w.attr = w.attr # allocation by individual attributes (random)
  return(P.vax)
}

def.params.net.city = function(N){
  # TODO: double check if / how seed should be used here
  P.net = list()
  P.net$N = N
  P.net$deg.shape = 0.255
  P.net$deg.rate = 0.032
  P.net$deg.shift = 0.913
  P.net$main.frac = .2
  P.net$main.w.deg.power = -1
  P.net$main.m.shape = 5
  P.net$main.m.rate = .2
  P.net$casu.m.shape = .3
  P.net$casu.m.rate = .3
  return(P.net)
}

make.net.city = function(P.net,city){
  # set.seed(P.net$seed) # TODO
  i = seqn(P.net$N)
  # sample degrees
  deg.i = round(rgamma(P.net$N,shape=P.net$deg.shape,P.net$deg.rate) + P.net$deg.shift)
  deg.i = degrees.balanced(deg.i)
  # generate main partners
  i.main = sample(i,round(P.net$main.frac*len(i)),p=deg.i^P.net$main.w.deg.power)
  ii.e.main = edges.random(i.main,shuffle=FALSE)
  r.e.main = round(rgamma(nrow(ii.e.main),shape=P.net$main.m.shape,rate=P.net$main.m.rate))
  ii.e.main.r = edges.repeated(ii.e.main,r.e.main)
  # generate casual partnerss
  deg.i.casu = deg.i
  deg.i.casu[i.main] = deg.i.casu[i.main] - 1
  ii.e.casu = edges.unloop(edges.from.degrees(i,deg.i.casu))
  r.e.casu = round(1+rgamma(nrow(ii.e.casu),shape=P.net$casu.m.shape,rate=P.net$casu.m.rate))
  ii.e.casu.r = edges.repeated(ii.e.casu,r.e.casu)
  # all contacts
  ii.e = rbind(ii.e.main.r,ii.e.casu.r)
  # attributes
  g.attr = list()
  g.attr$city = city
  i.attr = list()
  i.attr$city = rep(city,P.net$N)
  i.attr$par.p6m = deg.i
  e.attr = list()
  if (.debug){
    i.attr$sex.p6m = degrees.from.edges(i,ii.e)
    i.attr$main = factor(i %in% i.main,c(T,F),c('Yes','No'))
    e.attr$type = factor(c(rep('main',nrow(ii.e.main.r)),rep('casu',nrow(ii.e.casu.r))))
    e.attr$city = rep(city,nrow(ii.e))
    kr.e = index.repeated.edges(ii.e)
    e.attr$sex.index = kr.e[,1]
    e.attr$sex.total = kr.e[,2]
  }
  # graph object
  G = graph.obj(ii.e=ii.e,i=i,deg.i=deg.i,g.attr=g.attr,i.attr=i.attr,e.attr=e.attr)
  if (.debug){ G$attr$g$layout = graph.layout(G) }
  return(G)
}