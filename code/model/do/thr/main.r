source('utils/all.r')
source('graph/graph.r')
source('model/meta.r')
source('model/params.r')
source('model/epidemic.r')
source('model/plot.r')

.debug = FALSE

thr = list()
thr$tf = 180
# base
thr$base$N       = 10000
thr$base$N.s.epi = 10
thr$base$N.s.thr = 10
thr$base$p.I0    = .01
thr$base$vax.eff = .85
thr$base$vax.pri = 'random'

t.vec = epi.t(tf=thr$tf)

thr.fname.fig   = function(slug){root.path('out','fig','thr',slug,create=TRUE)}
thr.fname.rdata = function(slug){root.path('data','.rdata','thr',paste0(slug,'.rdata'),create=TRUE)}

thr.pri.wts = function(case,i.S,E){
  # define vaccine priority weights
  weights = switch(case,
    'random' = rep(1,len(i.S)),
    'degree' = E$P$G$attr$i$par[i.S])
}

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

thr.out = function(E,args){
  # estimate the vaccine thresholds (and Re) by generation for one model run
  # we do multiple estimates due to random distribution
  tree = thr.tree.gen(as.data.frame(E$tree)) # add generation column
  out = do.call(rbind,lapply(1:max(tree$gen),function(g){ # this generation:
    tree.g = tree[tree$gen==g,] # transmission pairs
    i.S = E$P$G$i[-tree$chi[tree$gen<g]]  # i of remaining susceptibles
    w.S = thr.pri.wts(args$vax.pri,i.S,E) # weights of remaining susceptibles
    n.S = length(i.S)                     # number of remaining susceptibles
    b.inf = i.S %in% tree.g$chi           # preventable infections
    n.par = length(unique(tree.g$par))    # number transmitting
    n.chi = length(tree.g$chi)            # number exposed
    n.thr = n.chi - n.par                 # infections we must avert
    thr = sapply(1:args$N.s.thr,function(s){  # random vaccine distribution & vax.eff
      b.inf.s = sample.wtd(b.inf & (runif(n.S) < args$vax.eff), w.S) # preventable infections
      thr.s = which(cumsum(b.inf.s) >= n.thr)[1] / n.S # find coverage needed for Re = 1
    })
    out.g = data.frame(gen=g,Re=n.chi/n.par,thr=thr) # collect outputs
  }))
}

thr.out.grid = function(args){
  # run the model for a grid of conditions, and compute thr.out, etc.
  args.grid = expand.grid(list.update(args, # clean agrs & make grid
    seed = 1:args$N.s.epi,
    N.I0 = round(args$N*args$p.I0)))
  out.grid = do.call(rbind,lapply(1:nrow(args.grid),function(i){ # for each combo:
    args.i = args.grid[i,]; print(args.i) # slice & print args
    E = epi.run(kw.call(def.params,args.i),t.vec) # run the epidemic
    out.i = cbind(args.i[c('vax.eff','vax.pri','seed')],thr.out(E,args.i)) # collect outputs
  }))
}

if (sys.nframe() == 0){
  args = thr$base
  out.grid = thr.out.grid(thr$base)
  out.grid$thr[is.na(out.grid$thr)] = 1
  g = ggplot(out.grid[out.grid$gen<=10,],aes(x=factor(gen),y=thr)) +
    geom_boxplot() +
    labs(x='Generation',y='Threshold') +
    theme_light()
  print(g)
}
