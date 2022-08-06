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
P.def = list(seed=5,N=50,N.I0=5,N.V10=0,N.V20=0)
G  = kw.call(def.params,P.def)$G
R0 = epi.run(kw.call(def.params,P.def,N.V10=0,G=G),t.vec)
R1 = epi.run(kw.call(def.params,P.def,N.V10=5,G=G),t.vec)
# update G$attr
G$attr$i$t = rep(t.vec,each=G$N.i)
G$attr$e$t = rep(t.vec,each=G$N.e)
G$attr$i$health.0 = c(t(R0$A[1:tf,]))
G$attr$i$health.1 = c(t(R1$A[1:tf,]))
M$health.0 = M$health; M$health.0$title = 'Base Case'
M$health.1 = M$health; M$health.1$title = 'With Vaccine'
# build & save gif
g = plot.network(G,list(fill='health.0',color='health.1',stroke=2,size=3)) +
  labs(title='Day {frame_time}') +
  transition_time(t) +
  ease_aes('linear')
animate(g,nframes=tf,fps=3,w=1200,h=1000,res=150,renderer=gifski_renderer('rfix-val.gif'))
