
def.params = function(seed=NULL,N=1000,...){
  set.seed(seed)
  P = list()
  # independent parameters (mostly)
  P$seed           = seed
  P$N              = N # pop size total
  P$N.I0           = 10 # number initially infected
  P$exp.deg.I0     = 0 # degree-exponent weight for initially infected
  P$dur.EI.rfun    = r.fun(rlnorm,meanlog=2.09,sdlog=0.46,rmin=3,rmax=21) # incubation period
  P$dur.IR.rfun    = r.fun(rgamma,shape=36,scale=0.58,rmin=14,rmax=28) # infectious period
  P$dur.IH.rfun    = r.fun(rgamma,shape=1.23,scale=4.05,rmin=2,rmax=20) # non-isolated period
  P$dur.IH.scale   = 1 # relative duration of non-isolated period
  P$p.asymp        = .15 # proportion of cases asymptomatic (never isolate)
  P$beta           = .90 # probability of transmission (per contact)
  P$vax.eff.dose   = c(.85,.88) # vaccine effectiveness by dose
  P$N.V0           = c(.00,.00) * P$N # total number initially vaccinated by dose
  P$exp.deg.V0     = 0 # degree-exponent weight for initially vaccinated
  P$p.detect.t     = interp.fun(c(0,30),c(0,.85),pad=TRUE) # probability of detection vs t
  P = list.update(P,...) # override any of the above
  P = def.params.net(P)
  # conditional parameters
  if (is.null(P$G)){ P$G = make.net(P) } # generate the sexual network
  P$beta.health = P$beta * # transmission prob by health state (susceptibility)
    c('S'=1,'E'=0,'I'=0,'H'=0,'R'=0,'V1'=1-P$vax.eff.dose[1],'V2'=1-P$vax.eff.dose[2])
  P$seed.state = .Random.seed # current state
  return(P)
}

def.params.s = function(seeds,...,.par=TRUE){
  # run def.params for a number (or vector) of seeds, parallel because net gen is expensive
  if (len(seeds)==1){ seeds = seqn(seeds) }
  P.s = par.lapply(seeds,def.params,...,.par=.par)
}

even = function(x){ x = round(x); x + x %% 2 }
rcap = function(x,xmin=1,xmax=Inf){ pmax(xmin,pmin(xmax,round(x))) }

def.params.net = function(P){
  P$t.vec = 1:180
  P$dur.main.rfun = r.fun(rweibull,shape=.521,scale=916,rmin=7,rmax=3650) # TODO
  P$dur.casu.rfun = r.fun(rweibull,shape=.5,  scale= 15,rmin=1,rmax=3650) # TODO
  P$dur.once.rfun = r.fun(rep,x=1)
  P$wt.rfun  = r.fun(runif)
  P$deg.once = 7   # TODO
  P$deg.casu = 2   # TODO
  P$f.sex    = 1/7 # TODO
  P$N.e.type = even(P$N / 2 * c(
    excl = .15/.607,
    open = .29/.607,
    casu = P$deg.casu,
    once = P$deg.once))
  return(P)
}

tt.fun = function(dur){
  t0 = round(runif(len(dur),-dur,180))
  tt = cbind(t0=t0,tf=t0+round(dur))
}

make.net = function(P){
  # the sexual network reflects all partnerships for 180 days
  i = seqn(P$N)
  w.i = P$wt.rfun(P$N)  # weight of forming a given partnership
  N.e = sum(P$N.e.type) # total partnerships
  # sample t0 & tf for all partnerships
  tt.excl = tt.fun(P$dur.main.rfun(P$N.e.type[1])) # excl
  tt.misc = rbind( # misc = open, casu, once
    tt.fun(P$dur.main.rfun(P$N.e.type[2])), # open
    tt.fun(P$dur.casu.rfun(P$N.e.type[3])), # casu
    tt.fun(P$dur.once.rfun(P$N.e.type[4]))) # once
  # assign excl
  i.excl  = sample(i,P$N.e.type[1]*2) # inds with 1 excl
  i.noex  = i[-i.excl]                # inds with no excl
  ii.excl = matrix(i.excl,ncol=2)     # excl pairs
  # assign misc pairs for compatible inds
  b.comp  = outer(tt.excl[,2],tt.misc[,1],`<`) | outer(tt.excl[,1],tt.misc[,2],`>`)
  ii.misc = t(sapply(1:nrow(tt.misc),function(e){ # for each misc partnership
    i = c(ii.excl[b.comp[,e],],i.noex) # available individuals
    sample(i,2,prob=w.i[i]) # sample individuals with weights
  }))
  # collect all partnerships
  ii.e = rbind(ii.excl, ii.misc)
  deg.i = tabulate(ii.e,P$N)
  # attributes
  g.attr = list()
  g.attr$dur = P$net.dur
  i.attr = list()
  i.attr$deg = deg.i
  e.attr = list()
  e.attr$t0  = c(tt.excl[,1],tt.misc[,1])
  e.attr$tf  = c(tt.excl[,2],tt.misc[,2])
  e.attr$dur = e.attr$tf - e.attr$t0
  if (.debug){ # expensive / not required
    # i.attr$stat = f.stat.i # TODO
    e.attr$type = factor(rep(names(P$N.e.type),P$N.e.type))
    # hist(deg.i,max(deg.i)) # DEBUG
  }
  # graph object
  G = graph.obj(ii.e=ii.e,i=i,deg.i=deg.i,g.attr=g.attr,i.attr=i.attr,e.attr=e.attr)
  # TODO: this results in .Random.seed depends on .debug: maybe move this after .Random.seed saved
  if (.debug & G$N.i < 1000){ G$attr$g$layout = graph.layout.fr(G) } # pre-compute consistent layout
  return(G)
}
