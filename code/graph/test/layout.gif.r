source('utils/all.r')
source('graph/graph.r')
library('gganimate')
# setup-graph
set.seed(0)
K = 50
N = 100
i = 1:N
deg.i = round(rexp(N,1) + 1)
ii.e = edges.unloop(edges.from.degrees(i,deg.i))
G = graph.obj(ii.e,i)
gif.pos = graph.layout.fr(G,niter=K,gif.pos=list())
G$attr$i$deg.i = factor(deg.i)
G$attr$i$k = rep(1:K,each=G$N.i)
G$attr$e$k = rep(1:K,each=G$N.e)
G$attr$g$layout = do.call(rbind,lapply(gif.pos,function(pos){
  apply(pos,2,scale)
}))
G$ii.e = do.call(rbind,lapply(1:K,function(k){
  G$ii.e + (k-1) * N
}))
g = plot.graph(G,list(fill='deg.i',size=3),list(alpha=.3)) +
  labs(title='{frame_time}') +
  transition_time(k) +
  ease_aes('linear') +
  scale_fill_viridis(discrete=TRUE,option='inferno')
animate(g,nframes=K,fps=3,w=800,h=800,res=150,renderer=gifski_renderer('layout.gif'))