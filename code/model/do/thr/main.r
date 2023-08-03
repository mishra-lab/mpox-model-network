source('utils/all.r')
source('graph/graph.r')
source('model/meta.r')
source('model/params.r')
source('model/epidemic.r')
source('model/plot.r')

.debug = FALSE

thr = list()
thr$tf = 180
thr$gen.max = 7
# base
thr$base$N       = 10000
thr$base$N.z.epi = 100
thr$base$N.z.vax = 10
thr$base$N.I0    = 10
thr$base$exp.deg.I0  = 0
thr$base$exp.deg.vax = 0
thr$base$vax.eff = .70

t.vec = epi.t(tf=thr$tf)

mean.le.1    = function(Re){ mean(Re <= 1) }
thr.fname    = function(slug){ root.path('out','fig','thr',slug,create=TRUE) }
thr.fig.save = function(slug,...){ fig.save(thr.fname(slug),...) }
thr.arg.save = function(slug,args){ cat(list.str(args),file=thr.fname(paste0(slug,'.args')),sep='\n') }

thr.tree.gen = function(tree,i.gen=0,gen=0){
  # add generation column to the tree
  if (gen==0){ tree = cbind(tree,gen=NA) }
  b.gen = tree$par %in% i.gen
  if (any(b.gen)){
    tree[b.gen,]$gen = gen
    tree = thr.tree.gen(tree,tree[b.gen,]$chi,gen+1)
  }
  return(tree)
}

thr.out = function(E,args,q.step=.01,deg.max=180){
  q.vec   = seq(0,1,q.step)
  deg.vec = seq(1,deg.max,1)
  # output templates
  out.gz.tmp = data.frame( # per generation (g) & stochastic vax (z)
    gen = NA,
    var = c('thr',rep('Re',len(q.vec))),
    strat = c(0,q.vec),
    value = NA)
  out.g.tmp = data.frame( # per generation (g) only
    gen = NA,
    var = c('prev',rep('deg',deg.max),rep('deg.inf',deg.max)),
    strat = c(0,deg.vec,deg.vec),
    value = NA)
  tree = thr.tree.gen(as.data.frame(E$tree)) # add generation column to E$tree
  # for each: generation
  out = rbind.lapply(1:max(tree$gen),function(gen){
    out.gz.tmp$gen = gen
    out.g.tmp$gen  = gen
    tree.g = tree[tree$gen==gen,] # transmission pairs
    i.S = E$P$G$i[-tree$chi[tree$gen<gen]]    # i of remaining susceptibles
    vax.w.S = E$P$G$deg.i[i.S]^args$exp.deg.vax # vax weights of susceptibles
    n.S = len(i.S)                  # number of remaining susceptible
    n.par = len(unique(tree.g$par)) # number at risk of transmit (parent) this gen
    n.chi = len(tree.g$chi)         # number at risk of acquire (child) this gen
    b.inf = i.S %in% tree.g$chi     # boolean at risk of acquire this gen among i.S
    # for each: stochastic vax (z)
    out.gz = rbind.lapply(1:args$N.z.vax,function(z){
      n.inf.z = cumsum(sample.wtd(b.inf & (runif(n.S) < args$vax.eff), vax.w.S))
      Re.z = (n.chi - n.inf.z) / n.par  # Re
      thr.z = which(Re.z <= 1)[1] / n.S # min coverage needed for Re <= 1
      out.gz.tmp$value = c(thr.z, quantile(Re.z,1-q.vec,names=FALSE,na.rm=TRUE))
      return(out.gz.tmp)
    })
    prev = (E$P$G$N.i - n.S) / E$P$G$N.i
    deg.inf = tabulate(E$P$G$deg.i[-i.S],deg.max)       # degrees of prior infected
    deg = deg.inf + tabulate(E$P$G$deg.i[+i.S],deg.max) # degrees overall
    out.g.tmp$value = c(prev, deg, deg.inf/deg)
    return(rbind(out.g.tmp,out.gz))
  })
}

thr.out.grid = function(args){
  # run the model for a grid of conditions, and compute thr.out, etc.
  args.grid = expand.grid(list.update(args,seed=1:args$N.z.epi)) # clean agrs & make grid
  out.grid = rbind.lapply(1:nrow(args.grid),function(i){ # for each combo:
    args.i = args.grid[i,]; print(args.i) # slice & print args
    E = epi.run(kw.call(def.params,args.i),t.vec) # run the epidemic
    out.i = cbind(args.i[c('N.I0','exp.deg.I0','exp.deg.vax','vax.eff','seed')],thr.out(E,args.i))
  },.par=TRUE)
}

plot.Re = function(out.grid,facet=NULL,clr=NULL,...){
  f = formula.fun('value','strat',facet,clr,...)
  out.grid.Re = aggregate(f,row.select(out.grid,var='Re',gen=1:thr$gen.max),mean.le.1)
  g = ggplot(out.grid.Re,aes(y=100*value,x=100*strat,!!!ensyms(color=clr,fill=clr,...))) +
    facet_grid(facet) +
    geom_ribbon(aes(ymin=0,ymax=100*value),alpha=.1,color=NA) +
    geom_line() +
    labs(y='Probability of Re < 1 (%)',x='Vaccine Coverage (%)')
  g = plot.clean(g)
}

plot.deg = function(out.grid,facet=NULL,clr=NULL,...){
  f = formula.fun('value','strat',facet,clr,...)
  out.grid.Re = aggregate(f,row.select(out.grid,var='deg.inf',gen=1:thr$gen.max),mean)
  g = ggplot(out.grid.Re,aes(y=100*value,x=strat,!!!ensyms(color=clr,fill=clr,...))) +
    facet_grid(facet) +
    geom_ribbon(aes(ymin=0,ymax=100*value),alpha=.1,color=NA) +
    geom_line() +
    labs(y='Proportion Infected (%)',x='P6M Partners')
  g = plot.clean(g)
}

plot.prev = function(out.grid,facet=NULL,clr=NULL,...){
  out.grid.prev = row.select(out.grid,var='prev',gen=1:thr$gen.max)
  g = ggplot(out.grid.prev,aes(y=100*value,!!!ensyms(color=clr,fill=clr,...))) +
    facet_grid(facet) +
    geom_violin(alpha=.1) +
    labs(y='Cumulative Infections (%)')
  g = plot.clean(g)
}

fig.gen = function(){
  args = thr$base
  out.grid = thr.out.grid(args)
  out.grid = as.factor.cols(out.grid,'gen')
  plot.Re  (out.grid,clr='gen');         thr.fig.save('gen.Re',  w=4,h=3)
  plot.deg (out.grid,clr='gen');         thr.fig.save('gen.deg', w=4,h=3)
  plot.prev(out.grid,clr='gen',x='gen'); thr.fig.save('gen.prev',w=4,h=3)
  thr.arg.save('gen',args)
}

fig.ve = function(){
  args = thr$base
  args = list.update(thr$base,vax.eff=c(.5,.6,.7,.8,.9))
  out.grid = thr.out.grid(args)
  out.grid = as.factor.cols(out.grid,c('gen','vax.eff'))
  plot.Re(out.grid,facet='gen',clr='vax.eff'); thr.fig.save('ve.Re',w=4,h=8)
  thr.arg.save('ve',args)
}

fig.exp = function(){
  args = thr$base
  args = list.update(thr$base,exp.deg.I0=c(0,1),exp.deg.vax=c(0,1))
  out.grid = thr.out.grid(args)
  out.grid = as.factor.cols(out.grid,c('gen','exp.deg.I0','exp.deg.vax'))
  plot.Re  (out.grid,facet='gen',clr='exp.deg.I0',lty='exp.deg.vax'); thr.fig.save('exp.Re',w=4,h=8)
  plot.deg (out.grid,facet='gen',clr='exp.deg.I0');                   thr.fig.save('exp.deg',w=4,h=8)
  plot.prev(out.grid,facet='gen',clr='exp.deg.I0',x='exp.deg.I0');    thr.fig.save('exp.prev',w=4,h=8)
  thr.arg.save('exp',args)
}

fig.gen()
# fig.ve()
fig.exp()

# args = list.update(thr$base,N.I0=c(3,30,300))
