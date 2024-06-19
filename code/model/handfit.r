# a collection of plots & tools to hand fit the epidemic

source('utils/all.r')
source('graph/graph.r')
source('model/meta.r')
source('model/params.r')
source('model/epidemic.r')
source('model/plot.r')

.debug = TRUE

fname = function(slug){
  # e.g. .../code/.tmp/handfit/2022-01-01/{slug}
  root.path('out','fig','.tmp','handfit',Sys.Date(),slug,create=TRUE)
}

handfit.run = function(seed=0,...){
  # run the model once
  P = def.params(seed=seed,...)
  E = epi.run(P)
}

handfit.plot.durs = function(P){
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

handfit.plot.epidemic = function(E){
  # main standard plots for one run
  # TODO: maybe add multiple model runs by default
  g = plot.epidemic(epi.output.melt(E$out,E$P),select=list(health=c('E','I','H'))) +
    labs(y='Individuals')
    fig.save(fname('epidemic'),w=8,h=4)
  g = plot.epidemic(epi.output.melt(E$out,E$P),select=list(var='inc')) +
    labs(y='Incidence')
    fig.save(fname('incidence'),w=8,h=4)
}

handfit.plot.doubling = function(E){
  g = plot.doubling(E,bw=7,lag=7,clip=60) +
    scale_y_continuous(breaks=seq(-60,60,10),trans='nsqrt') +
    lims(x=c(0,120))
    fig.save(fname('doubling'),w=6,h=4)
}

handfit.plot.G.distrs = function(E){
  # fill by main partner
  g = plot.G.distr(E$P,'i','deg',seq(0,40),fill='main.any') + geom_col() +
    labs(y='Individuals',x='Partners in past 6 months',fill='Main\nPartner\nP6M')
    fig.save(fname('i-deg-main-any'),w=8,h=4)
  g = plot.G.distr(E$P,'i','deg',seq(0,40),fill='main.now') + geom_col() +
    labs(y='Individuals',x='Partners in past 6 months',fill='Main\nPartner\nCurrent')
    fig.save(fname('i-deg-main-now'),w=8,h=4)
  # fill by health at tf (180)
  E$P$G$attr$i$health = E$P$G$attr$i$health.tf
  g = plot.G.distr(E$P,'i','deg',seq(0,40),fill='health') + geom_col() +
    labs(y='Individuals',x='Partners in past 6 months')
    fig.save(fname('i-deg-health'),w=8,h=4)
}

handfit.plot.tree = function(E){
  # common stuff
  add.cmaps = function(g,label,...){ g = g +
    scale_color_viridis(...,discrete=TRUE) +
    scale_fill_viridis(...,discrete=TRUE) +
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

handfit.network.gif = function(seed=0,N=100,tf=100,fps=5){
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

if (sys.nframe() == 0){
  .debug = TRUE
  E = handfit.run()
  handfit.plot.durs(E$P)
  handfit.plot.doubling(E)
  handfit.plot.epidemic(E)
  handfit.plot.G.distrs(E)
  handfit.plot.tree(E)
  handfit.network.gif()
}
