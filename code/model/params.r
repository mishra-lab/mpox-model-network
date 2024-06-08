
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
  P$t.vec = 1:P$dur.net
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
  P$dur.net = 180
  P$dur.main.rfun = r.fun(rweibull,shape=.521,scale=916)
  # P$dur.casu.rfun = r.fun(rweibull,shape=.315,scale=159) # TEMP
  P$dur.casu.rfun = r.fun(rweibull,shape=.5,scale=15) # TEMP
  P$deg.excl.rfun = r.fun(rgamma,  shape=.437,scale=12.7)
  P$deg.open.rfun = r.fun(rgamma,  shape=.461,scale=18.7)
  P$deg.noma.rfun = r.fun(rgamma,  shape=.527,scale=21.2)
  P$sex.rfun = function(d){ rcap(runif(len(d),0,pmin(d,30))) } # TEMP
  P$n.stat = even(P$N * c(excl=.15,open=.29,noma=NA) / .607)
  P$n.stat[3] = P$N - sum(P$n.stat[1:2])
  return(P)
}

make.net = function(P){
  # the sexual network reflects all partnerships for 180 days
  i = seqn(P$N)
  # stratify pop: if have excl, open, or noma (no main) during 180 days
  i.excl = 1:P$n.stat[1]
  i.open = 1:P$n.stat[2] + P$n.stat[1]
  i.noma = 1:P$n.stat[3] + sum(P$n.stat[1:2])
  f.stat.i = rep(names(P$n.stat),P$n.stat)
  # form main.excl
  ii.e.excl  = edges.random(i.excl)
  n.excl     = nrow(ii.e.excl)
  dur.e.excl = P$dur.main.rfun(n.excl)
  t0.e.excl  = runif(n.excl,-dur.e.excl,P$dur.net)
  tf.e.excl  = t0.e.excl + dur.e.excl
  # form main.open
  ii.e.open  = edges.random(i.open)
  n.open     = nrow(ii.e.open)
  dur.e.open = P$dur.main.rfun(n.open)
  t0.e.open  = runif(n.open,-dur.e.open,P$dur.net)
  tf.e.open  = t0.e.open + dur.e.open
  # count casu
  deg.casu.excl.i = rcap(P$deg.excl.rfun(P$n.stat[1])-1,0) # unused
  deg.casu.open.i = rcap(P$deg.open.rfun(P$n.stat[2])-1,0)
  deg.casu.noma.i = rcap(P$deg.noma.rfun(P$n.stat[3]))
  n.casu = c(sum(deg.casu.excl.i),sum(deg.casu.open.i),sum(deg.casu.noma.i))
  i.casu.noex = sample(c( # shuffled indices of noex (non-excl)
    rep(i.open,deg.casu.open.i),
    rep(i.noma,deg.casu.noma.i)))
  # define gaps for casu.excl (ce) around main.excl
  e = order(ii.e.excl) %% n.excl; e[e==0] = n.excl # order of *.e.excl for i.excl
  dt.pre = pmax(0,          t0.e.excl[e]) # gap before excl
  dt.pos = pmax(0,P$dur.net-tf.e.excl[e]) # gap after excl
  ct.pre = cumsum(dt.pre) # pre-compute cumsum
  ct.pos = cumsum(dt.pos) # pre-compute cumsum
  dur.ce = P$dur.casu.rfun(n.casu[1]) # durations
  # sample dummy var (ctg) = cumulative gap time around main.excl
  ctg.ce = runif(n.casu[1],-sum(dt.pre),+sum(dt.pos))
  b.pre = ctg.ce < 0 # casu before excl
  ctg.ce.pre = -ctg.ce[ b.pre] # cumulative tf before excl
  ctg.ce.pos = +ctg.ce[!b.pre] # cumulative t0 after excl
  # find which excl each casu is before/after:
  i.ce.pre = findInterval(ctg.ce.pre,ct.pre)+1
  i.ce.pos = findInterval(ctg.ce.pos,ct.pos)+1
  # define valid t0 & tf from the above
  tf.ce.pre =              ct.pre[i.ce.pre] - ctg.ce.pre
  t0.ce.pos = P$dur.net - (ct.pos[i.ce.pos] - ctg.ce.pos)
  t0.ce.pre = tf.ce.pre - dur.ce[ b.pre]
  tf.ce.pos = t0.ce.pos + dur.ce[!b.pre]
  i1.e.ce = c( i.ce.pre, i.ce.pos)
  t0.e.ce = c(t0.ce.pre,t0.ce.pos)
  tf.e.ce = c(tf.ce.pre,tf.ce.pos)
  # form casu among excl-noex
  i.e.ce = c(i1.e.ce, i.casu.noex[1:n.casu[1]])
  ii.e.ce = edges.random(i.e.ce,shuffle=FALSE)
  deg.casu.excl.i = tabulate(i.e.ce,P$n.stat[1])
  # form casu among noex-noex
  ii.e.cx  = edges.unloop(edges.random(i.casu.noex[-(1:n.casu[1])]))
  n.cx     = nrow(ii.e.cx)
  dur.e.cx = P$dur.casu.rfun(n.cx)
  t0.e.cx  = runif(n.cx,-dur.e.cx,P$dur.net)
  tf.e.cx  = t0.e.cx + dur.e.cx
  # all ptrs
  ii.e = rbind(ii.e.excl, ii.e.open, ii.e.ce, ii.e.cx)
  t0.e = rcap(c(t0.e.excl, t0.e.open, t0.e.ce, t0.e.cx),1,180)
  tf.e = rcap(c(tf.e.excl, tf.e.open, tf.e.ce, tf.e.cx),1,180)
  dur.e = tf.e - t0.e
  sex.e = P$sex.rfun(dur.e)
  deg.i = c(1+deg.casu.excl.i,1+deg.casu.open.i,deg.casu.noma.i)
  # attributes
  g.attr = list()
  g.attr$dur = P$net.dur
  i.attr = list()
  i.attr$deg = deg.i
  e.attr = list()
  e.attr$t0  = t0.e
  e.attr$tf  = tf.e
  e.attr$dur = dur.e
  e.attr$sex = sex.e
  if (.debug){ # expensive / not required
    i.attr$sex = aggregate(sex~i,cbind(i=c(ii.e),sex=sex.e),sum)$sex
    i.attr$stat = f.stat.i
    e.attr$type = factor(c(
      rep('excl',nrow(ii.e.excl)),
      rep('open',nrow(ii.e.open)),
      rep('casu',nrow(ii.e.ce)+nrow(ii.e.cx))))
  }
  # graph object
  G = graph.obj(ii.e=ii.e,i=i,deg.i=deg.i,g.attr=g.attr,i.attr=i.attr,e.attr=e.attr)
  # TODO: this results in .Random.seed depends on .debug: maybe move this after .Random.seed saved
  if (.debug & G$N.i < 1000){ G$attr$g$layout = graph.layout.fr(G) } # pre-compute consistent layout
  return(G)
}
