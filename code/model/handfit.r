# a collection of plots & tools to hand fit the epidemic

source('utils/all.r')
source('graph/graph.r')
source('model/meta.r')
source('model/params.r')
source('model/epidemic.r')
source('model/plot.r')

.debug = TRUE
t = epi.t(tf=180)

fname = function(slug){
  # e.g. .../code/.tmp/handfit/2022-01-01/{slug}
  root.path('code','.tmp','handfit',Sys.Date(),slug,create=TRUE)
}

handfit.run = function(seed=0,...){
  # run the model once
  P = def.params(seed=seed,...)
  R = epi.results(P,t,epi.run(P,t))
}

handfit.plot.epidemic = function(R){
  # main standard plots for one run
  # TODO: maybe add multiple model runs by default
  plot.epidemic(epi.output.melt(R$out,R$P)) +
    labs(y='Individuals')
    fig.save(fname('epidemic'),w=8,h=4)
  plot.epidemic(epi.output.melt(R$out,R$P),select=list(var='inc')) +
    labs(y='Incidence')
    fig.save(fname('incidence'),w=8,h=4)
}

handfit.plot.distribs = function(R){
  # fill by main partner
  plot.distribution(R$P,'i','par',seq(1,20),fill='main') + geom_col() +
    labs(y='Individuals',x='Partners in past 6 months',fill='Main\nPartner\nP6M')
    fig.save(fname('i-par-main'),w=8,h=4)
  plot.distribution(R$P,'i','sex',seq(0,100,10),fill='main') + geom_col() +
    labs(y='Individuals',x='Sex in past 6 months',fill='Main\nPartner\nP6M')
    fig.save(fname('i-sex-main'),w=8,h=4)
  # fill by health at tf (180)
  plot.distribution(R$P,'i','par',seq(1,20),fill='health') + geom_col() +
    labs(y='Individuals',x='Partners in past 6 months')
    fig.save(fname('i-par-health'),w=8,h=4)
  plot.distribution(R$P,'i','sex',seq(0,100,10),fill='health') + geom_col() +
    labs(y='Individuals',x='Sex in past 6 months')
    fig.save(fname('i-sex-health'),w=8,h=4)
}

handfit.network.gif = function(seed=0,N=100,tf=100,fps=10){
  # make a gif from plot.network, showing health states over time
  library('gganimate')
  # define params & run model
  P = def.params(seed=seed,N=N,N.city=c(N,0,0),N.I0.city=c(5,0,0),vax.params.phase=list(),
    net.params.city=list(A=def.params.net.city(N)))
  t.vec = epi.t(tf=tf)
  A = epi.run(P,t.vec)
  # hijack attr: write too many values that will be recycled by plot.network
  P$G$attr$i$t = rep(t.vec,each=P$G$N.i)
  P$G$attr$e$t = rep(t.vec,each=P$G$N.e)
  P$G$attr$i$health = c(t(A[1:tf,]))
  # build & save gif
  g = plot.network(P$G,list(fill='health')) + labs(title='Day {frame_time}') +
    transition_time(t) +
    ease_aes('linear')
  f = fname(sprintf('gif-n%d-s%d-t%d.gif',N,seed,tf)) # filename
  animate(g,nframes=tf,fps=fps,w=600,h=500,renderer=gifski_renderer(f))
}

R = handfit.run()
handfit.plot.epidemic(R)
handfit.plot.distribs(R)
handfit.network.gif(seed=0,N=100,tf=100)
handfit.network.gif(seed=0,N=300,tf=100)





