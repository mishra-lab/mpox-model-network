# a collection of plots & tools to hand fit the epidemic

source('utils/all.r')
source('graph/graph.r')
source('model/meta.r')
source('model/params.r')
source('model/epidemic.r')
source('model/plot.r')
source('model/data.r')

fname = function(slug){
  # e.g. .../code/.tmp/fit/2022-01-01/{slug}
  root.path('out','fig','.tmp','fit',Sys.Date(),slug,create=TRUE)
}

# -----------------------------------------------------------------------------
# num ptrs in 6 months vs data

fit.p6m = function(N.s=21){
  X.dat = row.select(main.p6m(main.stat()),city='avg')
  fit = cbind(#    init,   min,  max
    w.shape    = c( 0.764, 0.2  , 9. ),
    r.ptr.casu = c( 0.015, 0.001, 0.3),
    r.ptr.once = c( 0.016, 0.001, 0.3),
    w.pwr.excl = c(+0.103,-1    ,+1  ),
    w.pwr.open = c(-0.115,-1    ,+1  ))
  err.fun = function(fit){
    X.mod = fit.p6m.mod(N.s=N.s,fit=as.list(fit))
    e = sum(abs(X.dat$cp - X.mod$cp)^2)
    print(c(fit,err=e))
    return(e) }
  opt = optim(unlist(fit[1,]),err.fun,method='L-BFGS-B',
    lower=unlist(fit[2,]), upper=unlist(fit[3,]))
  # opt = list(par=unlist(fit[1,])) # DEBUG
  print(opt)
  X.mod = fit.p6m.mod(N.s=1000,fit=as.list(opt$par))
  fit.p6m.plot(X.dat,X.mod)
}

fit.p6m.mod = function(N.s,...,x=0:300,N=1000){
  # get num ptrs in 6 months from N.s models
  P.s = def.params.s(N.s,N=N,t.max=180,...,.par=FALSE)
  X.mod = rbind.lapply(P.s,function(P){
    n = split(P$G$attr$i$n.ptr.tot,P$G$attr$i$main.any)
    X.mod.s = rbind.lapply(names(n),function(stat){
      nt = tabulate0(n[[stat]],max(x))
      data.frame(seed=P$seed,status=stat,x=x,p=nt/sum(nt),cp=cumsum(nt)/sum(nt))
    })
  })
}

fit.p6m.plot = function(X.dat,X.mod,x.max=50){
  facet = c(p='Proportion',cp='Cumulative Proportion')
  X = rbind(cbind(X.dat,seed=0,src='Data'),cbind(X.mod,city='mod',src='Model'))
  X = melt(X,measure=names(facet),var='facet')
  X$facet = factor(X$facet,names(facet),facet)
  g = ggplot(X,aes(x=x,y=value,color=status,fill=status,linetype=src)) +
    facet_grid('facet',scales='free_y') +
    stat_summary(geom='ribbon',fun.min=min,fun.max=max,alpha=.3,color=NA) +
    stat_summary(geom='line',fun=median) +
    lims(x=c(0,x.max)) +
    labs(y='Proportion',x='Total partners in 6 months',linetype='')
  g = add.meta.scales(g,list(color='status',fill='status'))
  g = plot.clean(g)
  fig.save(fname('p6m'),w=5,h=5)
}

# -----------------------------------------------------------------------------
# other param plots

fit.plot.durs = function(){
  P = def.params()
  funs = list(
    dur.EI = 'Incubation Period',
    dur.IR = 'Infectious Period',
    dur.IH = 'Non-Isolated Period')
  clrs = unname(M$health$color[2:4])
  rfuns = setNames(P[paste0(names(funs),'.rfun')],funs)
  g = plot.durs(rfuns) +
    scale_color_manual(values=clrs) +
    scale_fill_manual(values=clrs)
    fig.save(fname('rfuns'),w=8,h=5)
}

# -----------------------------------------------------------------------------
# model outputs plots

fit.run.epi = function(seed=0,...){
  # run the model once
  P = def.params(seed=seed,...)
  E = epi.run(P)
}

fit.plot.epidemic = function(E){
  # main standard plots for one run
  # TODO: maybe add multiple model runs by default
  g = plot.epidemic(epi.output.melt(E),select=list(health=c('E','I','H'))) +
    labs(y='Individuals')
    fig.save(fname('epidemic'),w=8,h=4)
  g = plot.epidemic(epi.output.melt(E),select=list(var='inc')) +
    labs(y='Incidence')
    fig.save(fname('incidence'),w=8,h=4)
}

fit.plot.tree = function(E){
  # common stuff
  add.cmaps = function(g,label,...){ g = g +
    viridis::scale_color_viridis(...,discrete=TRUE) +
    viridis::scale_fill_viridis(...,discrete=TRUE) +
    labs(color=label,fill=label) }
  # transmission by generation
  E$tree$gen.f = factor(E$tree$gen)
  g = plot.tree(E$tree,fill='gen.f')
    g = add.tree.margin(add.cmaps(g,'Generation',option='inferno',begin=.2,end=.9))
    fig.save(fname('tree-gen'),g=g,w=8,h=6)
  # transmission by n children
  E$tree$n.chi.cut = int.cut(E$tree$n.chi.dir,c(0,1,2,6))
  g = plot.tree(E$tree,fill='n.chi.cut')
    g = add.tree.margin(add.cmaps(g,'Transmission',option='viridis'))
    fig.save(fname('tree-n-child'),g=g,w=8,h=6)
}

# -----------------------------------------------------------------------------
# dynamic network gif (expensive)

fit.network.gif = function(seed=0,N=100,tf=100,fps=5){
  # make a gif from plot.network, showing health states over time
  library('gganimate')
    # define params & run model
  P = def.params(seed=seed,N=N,N.I0=5)
  E = epi.run(P)
  t.vec = P$t.vec[1:tf]
  # hijack attr: write too many values that will be recycled by plot.network
  act.e = outer(t.vec,P$G$attr$e$t0,`>=`) & outer(t.vec,P$G$attr$e$tf,`<`)
  P$G$attr$i$t = rep(t.vec,each=P$G$N.i)
  P$G$attr$e$t = rep(t.vec,each=P$G$N.e)
  P$G$attr$e$active = factor(c(t(act.e)))
  P$G$attr$i$health = c(t(E$A[1:tf,]))
  # build & save gif
  g = plot.network(P$G,list(fill='health'),list(linetype='active',alpha=1)) +
    labs(title='Day {frame_time}') +
    transition_time(t) +
    ease_aes('linear')
  f = fname(sprintf('gif-n%d-s%d-t%d.gif',N,seed,tf)) # filename
  animate(g,nframes=tf,fps=fps,w=600,h=500,renderer=gifski_renderer(f))
}

# -----------------------------------------------------------------------------
# main

if (sys.nframe() == 0){
  .debug = FALSE
  fit.p6m()
  fit.plot.durs()
  .debug = TRUE
  E = fit.run.epi()
  fit.plot.epidemic(E)
  fit.plot.tree(E)
  fit.network.gif()
}
