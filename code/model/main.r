source('utils/all.r')
source('graph/graph.r')
source('model/meta.r')
source('model/params.r')
source('model/epidemic.r')
source('model/plot.r')
library('ggplot2')

.debug = FALSE

t = epi.t(tf=180)

# DEBUG: run one
# P = def.params(seed=0)
# R = epi.results(P,t,epi.run(P,t))
# g = plot.epidemic(epi.output.melt(R$out,P),select=list(var='inc',health=c('S','V1'))); print(g)
# plot.network(R$P$G,list(fill='health'),list()); fig.save('.tmp/network',w=7,h=7)

# build + run many
N.s = 8
P.s = def.params.s(N.s)
R.s = epi.run.s(P.s,t)
out.long.s = epi.output.melt.s(R.s,P.s)

# plot prevalence by city
g = plot.epidemic(out.long.s,select=list(var='prev',city=c('A'))) +
  labs(y='Prevalence',color='State',fill='State')
fig.save('.tmp/prev.city',w=8,h=4)

# plot networks
par.lapply(seq(N.s),function(s){
  plot.network(R.s[[s]]$P$G,list(fill='health',color='inf.src',shape='inf.src'))
  fig.save('.tmp/net-',s,w=7,h=7)
})
