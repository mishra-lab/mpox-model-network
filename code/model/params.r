
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

rfun.sex = function(n,q.dur,q.freq,cor=-.9,dmax=180,fmax=.5){
  # draw and multiply neg-correlated samples from q.dur & q.freq
  # after bounding by dmax & fmax
  u    = pnorm(mvtnorm::rmvnorm(n,sigma=matrix(c(1,cor,cor,1),2)))
  dur  = pmin(dmax,q.dur [ceiling(u[,1]*length(q.dur))])
  freq = pmin(fmax,q.freq[ceiling(u[,2]*length(q.freq))])
  n.sex = pmax(1,round(dur * freq))
}

make.net = function(P){
  # the sexual network reflects all contacts (sex) occuring in 180 days
  # including multiple contacts per partnership
  i = seqn(P$N)
  # stratify pop & fix odd numbers
  n.stat = round(P$N*d.stat)
  n.odd  = n.stat[2:3] %% 2 # open,excl
  n.stat = n.stat + c(-sum(n.odd),n.odd)
  f.stat.i = rep(names(d.stat),n.stat)
  i.stat = split(i,f.stat.i)
  # get degrees
  deg.s = sapply(names(d.stat),function(s){
    pmax(1,sample(x=d.ptrs$n,size=n.stat[[s]],prob=d.ptrs[[s]],replace=TRUE)) })
  deg.i      = c(deg.s$casu, deg.s$open,   deg.s$excl)
  deg.i.casu = c(deg.s$casu, deg.s$open-1, deg.s$excl-1)
  # generate main ptrs
  ii.e.main = rbind(
    edges.random(i.stat$excl,shuffle=TRUE),
    edges.random(i.stat$open,shuffle=TRUE))
  sex.e.main = rfun.sex(nrow(ii.e.main),q.dur$main,q.freq$main)
  # generate casu ptrs
  ii.e.casu  = edges.unloop(edges.from.degrees(i,deg.i.casu))
  sex.e.casu = rfun.sex(nrow(ii.e.casu),q.dur$casu,q.freq$casu)
  # all ptrs
  ii.e = rbind(ii.e.main,ii.e.casu)
  sex.e = c(sex.e.main,sex.e.casu)
  # attributes
  g.attr = list()
  g.attr$dur = 180
  i.attr = list()
  i.attr$deg = deg.i
  e.attr = list()
  e.attr$sex = sex.e
  if (.debug){ # expensive / not required
    i.attr$sex = aggregate(sex~i,cbind(i=c(ii.e),sex=sex.e),sum)$sex
    i.attr$stat = f.stat.i
    e.attr$type = factor(c(rep('main',nrow(ii.e.main)),rep('casu',nrow(ii.e.casu))))
  }
  # graph object
  G = graph.obj(ii.e=ii.e,i=i,deg.i=deg.i,g.attr=g.attr,i.attr=i.attr,e.attr=e.attr)
  # TODO: this results in .Random.seed depends on .debug: maybe move this after .Random.seed saved
  if (.debug & G$N.i < 1000){ G$attr$g$layout = graph.layout.fr(G) } # pre-compute consistent layout
  return(G)
}

# pre-load data
priv.csv = function(fname){ read.csv(root.path('data','.private',paste0(fname,'.csv'))) }
d.ptrs = priv.csv('ptrs_distr')
d.stat = priv.csv('stat_distr')
q.freq = priv.csv('freq_quant')
q.dur  = priv.csv('dur_quant')
