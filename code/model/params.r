
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
  P$beta           = .67 # probability of transmission (per contact)
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
  P$w.ptr.rfun  = function(n){ runif(n) }
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

sample.tt = function(dur){
  # sample t0,tf given dur which guarantees all partnerships
  # will overlap with [1,180], i.e. t0 <= 180 or tf >= 1
  t0 = round(runif(len(dur),1-dur,180))
  tt = cbind(t0=t0,tf=t0+round(dur))
}

e.match = function(ii.excl,ii.misc){
  # for each i in ii.misc, return the row where i is found in ii.excl, or NA
  # len(e) = 2*nrow(ii.misc) corresponding to c(ii[,1],ii[,2])
  e = match(ii.misc,ii.excl) %% nrow(ii.excl)
  e[e==0] = nrow(ii.excl)
  return(e)
}

assign.ii = function(tt.excl,tt.misc,ii.excl,i,w.i){
  # assign misc pairs (ii) among all i while avoiding overlaps with tt.excl
  if (nrow(tt.misc)==0){ return(NULL) } # end of recursion
  ii.misc = matrix(sample(i,nrow(tt.misc)*2,rep=TRUE,prob=w.i),ncol=2) # random pairs
  e = e.match(ii.excl,ii.misc) # find rows of ii.misc in ii.excl
  b.val.ee = matrix(ncol=2, # boolean: valid pairs for tt.excl vs tt.misc
    is.na(e) | tt.misc[,2] < tt.excl[e,1] | tt.misc[,1] > tt.excl[e,2])
  b.val.e = ii.misc[,1] != ii.misc[,2] & b.val.ee[,1] & b.val.ee[,2] # combine for ii[,1] & ii[,2]
  # recursively re-attempt invalid pairs
  ii.misc[!b.val.e,] = assign.ii(tt.excl,tt.misc[!b.val.e,,drop=FALSE],ii.excl,i,w.i)
  return(ii.misc)
}

make.net = function(P){
  # the sexual network reflects all partnerships for 180 days
  i = seqn(P$N)
  w.i = P$w.ptr.rfun(P$N) # weight for forming non-excl partnerships
  N.e = sum(P$N.e.type) # total partnerships
  # sample t0 & tf for all partnerships
  tt.excl = sample.tt(P$dur.main.rfun(P$N.e.type[1]))
  tt.open = sample.tt(P$dur.main.rfun(P$N.e.type[2]))
  tt.casu = sample.tt(P$dur.casu.rfun(P$N.e.type[3]))
  tt.once = sample.tt(P$dur.once.rfun(P$N.e.type[4]))
  # assign pairs
  i.excl = sample(i,P$N.e.type[1]*2,prob=1-w.i) # inds with 1 excl-main
  i.open = sample(i[-i.excl],P$N.e.type[2]*2)   # inds with 1 open-main
  ii.excl = matrix(i.excl,ncol=2) # excl pairs
  ii.open = matrix(i.open,ncol=2) # open pairs
  ii.casu = assign.ii(tt.excl,tt.casu,ii.excl,i,w.i) # casu pairs
  ii.once = assign.ii(tt.excl,tt.once,ii.excl,i,w.i) # once pairs
  ii.e = rbind(ii.excl, ii.open, ii.casu, ii.once) # all pairs
  # attributes
  g.attr = list()
  g.attr$dur = P$net.dur
  i.attr = list()
  i.attr$deg = tabulate(ii.e,P$N)
  e.attr = list()
  e.attr$t0  = c(tt.excl[,1],tt.open[,1],tt.casu[,1],tt.once[,1])
  e.attr$tf  = c(tt.excl[,2],tt.open[,2],tt.casu[,2],tt.once[,2])
  e.attr$dur = e.attr$tf - e.attr$t0
  if (.debug){ # expensive / not required
    i.attr$w.ptr = w.i
    i.attr$stat = as.factor(ifelse(i %in% ii.excl,'excl',
                            ifelse(i %in% ii.open,'open','noma')))
    e.attr$type = factor(rep(names(P$N.e.type),P$N.e.type))
    # hist(i.attr$deg,max(i.attr$deg)) # DEBUG
  }
  # graph object
  G = graph.obj(ii.e=ii.e,i=i,g.attr=g.attr,i.attr=i.attr,e.attr=e.attr)
  # TODO: this results in .Random.seed depends on .debug: maybe move this after .Random.seed saved
  if (.debug & G$N.i < 1000){ G$attr$g$layout = graph.layout.fr(G) } # pre-compute consistent layout
  return(G)
}
