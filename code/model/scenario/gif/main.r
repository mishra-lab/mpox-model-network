source('utils/all.r')
source('graph/graph.r')
source('model/meta.r')
source('model/params.r')
source('model/epidemic.r')
source('model/plot.r')
library('gganimate')
# define params & run models
.debug = TRUE
tf = 60
t.vec = epi.t(tf=tf)
P.def = list(seed=0,N=100,N.I0=3,N.V0=c(0,0)) # default params
G  = kw.call(def.params,P.def)$G # pre-compute fixed network
E = list( # run model with & without vax
  'base' = epi.run(kw.call(def.params,P.def,N.V0=c( 0,0),G=G),t.vec),
  'vax'  = epi.run(kw.call(def.params,P.def,N.V0=c(25,0),G=G),t.vec))
# replicate the tree data for each scenario & day, and copy health states
tree = do.call(rbind,lapply(names(E),function(scen){ # scenarios
  tree = tree.pc.map(E$base$tree,c('t','pos')) # for geom_segment
  tree.t = do.call(rbind,lapply(t.vec,function(tj){ # days (tj)
    cbind(tree,'health'=c(NA,E[[scen]]$A[tj,tree$chi]),'t.vec'=tj,'scen'=scen)
  }))
}))
# clean-up some data
tree$scen = factor(tree$scen,labels=c('No Vaccine','With Vaccine')) # prettify
t.data = data.frame(t.vec=t.vec,p0=-0.1,p1=1.1) # dummy hline data (geom_hline had issues)
# plot & generate gif
g = plot.tree(tree,size=2,fill='health',pc.map=FALSE) +
  geom_segment(data=t.data,aes(y=t.vec,yend=t.vec,x=p0,xend=p1),alpha=.05,size=1) +
  facet_wrap('scen') +
  labs(title='Day {frame_time}') +
  transition_time(t.vec) +
  ease_aes('linear')
g = add.meta.scales(g,list(color='health',fill='health'))
animate(g,nframes=tf,fps=3,w=800,h=500,res=150,renderer=gifski_renderer('rfix-val.gif'))

