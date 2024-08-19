
def.params = function(seed=NULL,N=1000,context='engage',t.max=180,...){
  set.seed(seed)
  P = list()
  # independent parameters (mostly)
  P$seed           = seed
  P$context        = as.character(context)
  P$t.max          = t.max
  P$N              = N # pop size total
  P$N.I0           = 10 # number initially infected
  P$I0.pfun        = function(P){ P$G$attr$i$w.ptr } # probs for initially infected
  P$dur.EI.rfun    = r.fun(rlnorm,meanlog=2.09,sdlog=0.46,rmin=3,rmax=21) # incubation period
  P$dur.IR.rfun    = r.fun(rgamma,shape=36,scale=0.58,rmin=14,rmax=28) # infectious period
  P$dur.IH.rfun    = r.fun(rgamma,shape=1.23,scale=4.05,rmin=2,rmax=20) # non-isolated period
  P$dur.IH.scale   = 1 # relative duration of non-isolated period
  P$p.asymp        = .15 # proportion of cases asymptomatic (never isolate)
  P$beta           = .67 # probability of transmission (per contact)
  P$vax.eff.dose   = c(.85,.88) # vaccine effectiveness by dose
  P$N.V0           = c(.00,.00) * P$N # total number initially vaccinated by dose
  P$V0.pfun        = function(P){ rep(1,P$N) } # probs of initially vaccinated
  P$fit = def.params.fit(P$context)
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

def.params.fit = function(context){
  fit = list(
  'engage' = list(
    w.sdlog      =  1.3,   # fit
    r.ptr.casu   =  0.020, # fit
    r.ptr.once   =  0.020, # fit
    w.pwr.excl   = +0.200, # fit
    w.pwr.open   = -0.300, # fit
    d.main.scale =   1059, # data
    d.casu.scale =    299, # data
    p.excl       =  0.154, # data
    p.open       =  0.309),# data
  'kenya' = list(
    w.sdlog      = 0.4,    # fit
    r.ptr.casu   = 0.005,  # fit
    r.ptr.once   = 0.030,  # fit
    w.pwr.excl   = 0,      # fit
    w.pwr.open   = 0,      # fit
    d.main.scale = 720,    # fit
    d.casu.scale =  60,    # fit
    p.excl       = 0.1,    # fit
    p.open       = 0.5)    # fit
  )[[context]]
}

def.params.net = function(P){
  P$t.vec         = 1:P$t.max
  P$dur.main.rfun = r.fun(rweibull,shape=.585,scale=P$fit$d.main.scale,rmin=1,rmax=3650)
  P$dur.casu.rfun = r.fun(rweibull,shape=.387,scale=P$fit$d.casu.scale,rmin=1,rmax=3650)
  P$dur.once.rfun = r.fun(rep,x=1)
  P$w.ptr.rfun  = function(n){ rlnorm(n,meanlog=1,sdlog=P$fit$w.sdlog) }
  P$f.sex   = .25 # TODO
  p.main.xs = 1-.046*P$t.max^.39 # handfit
  P$N.e.type = even(P$N / 2 * c(
    excl = P$fit$p.excl/p.main.xs,
    open = P$fit$p.open/p.main.xs,
    casu = P$fit$r.ptr.casu*P$t.max,
    once = P$fit$r.ptr.once*P$t.max))
  return(P)
}

sample.tt = function(durs,t.max){
  # sample t0,tf given dur which guarantees all partnerships
  # will overlap with [1,t.max], i.e. t0 <= t.max or tf >= 1
  t0 = round(runif(len(durs),1-durs,t.max))
  tt = cbind(t0=t0,tf=t0+round(durs))
}

e.match = function(ii.excl,ii.misc){
  # for each i in ii.misc, return the row where i is found in ii.excl, or NA
  # len(e) = 2*nrow(ii.misc) corresponding to c(ii) = c(ii[,1],ii[,2])
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
  # the sexual network reflects all partnerships for P$t.max days
  i = seqn(P$N)
  w.i = P$w.ptr.rfun(P$N) # weight for forming non-excl partnerships
  N.e = sum(P$N.e.type) # total partnerships
  # sample t0 & tf for all partnerships
  tt.excl = sample.tt(P$dur.main.rfun(P$N.e.type[1]),P$t.max)
  tt.open = sample.tt(P$dur.main.rfun(P$N.e.type[2]),P$t.max)
  tt.casu = sample.tt(P$dur.casu.rfun(P$N.e.type[3]),P$t.max)
  tt.once = sample.tt(P$dur.once.rfun(P$N.e.type[4]),P$t.max)
  # assign pairs
  w.excl = w.i^P$fit$w.pwr.excl
  w.open = w.i^P$fit$w.pwr.open
  i.excl = sample.wtd(i,w=w.excl,n=P$N.e.type[1]*2) # inds with 1 excl-main
  i.open = sample.wtd(i[-i.excl],w=w.open[-i.excl],n=P$N.e.type[2]*2) # inds with 1 open-main
  ii.excl = matrix(i.excl,ncol=2) # excl pairs
  ii.open = matrix(i.open,ncol=2) # open pairs
  ii.casu = assign.ii(tt.excl,tt.casu,ii.excl,i,w.i) # casu pairs
  ii.once = assign.ii(tt.excl,tt.once,ii.excl,i,w.i) # once pairs
  ii.e = rbind(ii.excl, ii.open, ii.casu, ii.once) # all pairs
  # attributes
  g.attr = list()
  g.attr$dur = P$t.max
  e.attr = list()
  e.attr$t0  = c(tt.excl[,1],tt.open[,1],tt.casu[,1],tt.once[,1])
  e.attr$tf  = c(tt.excl[,2],tt.open[,2],tt.casu[,2],tt.once[,2])
  e.attr$dur = e.attr$tf - e.attr$t0
  e.attr$type = factor(rep(names(P$N.e.type),P$N.e.type))
  i.attr = list()
  i.attr$n.ptr.tot = tabulate(ii.e,P$N)
  i.attr$n.ptr.now = tabulate(ii.e[e.attr$t0 <= 1,],P$N)
  i.attr$w.ptr = w.i
  i.attr$main.any = factor(levels=M$main$name,
    ifelse(i %in% ii.excl,'excl',
    ifelse(i %in% ii.open,'open','noma')))
  i.attr$main.now = factor(levels=M$main$name,
    ifelse(i %in% ii.excl[tt.excl[,2] >= P$t.max,],'excl',
    ifelse(i %in% ii.open[tt.open[,2] >= P$t.max,],'open','noma')))
  # hist(i.attr$deg,max(i.attr$deg)) # DEBUG
  # graph object
  G = graph.obj(ii.e=ii.e,i=i,g.attr=g.attr,i.attr=i.attr,e.attr=e.attr)
  # TODO: this results in .Random.seed depends on .debug: maybe move this after .Random.seed saved
  if (.debug & G$N.i < 1000){ G$attr$g$layout = graph.layout.fr(G) } # pre-compute consistent layout
  return(G)
}
