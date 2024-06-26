source('utils/all.r')
source('graph/graph.r')
source('model/meta.r')
source('model/params.r')
source('model/epidemic.r')
source('model/plot.r')

.debug = TRUE

# DEBUG: run one
# P = def.params(seed=0)
# E = epi.run(P)
# plot.epidemic(epi.output.melt(E$out,P)); fig.save('.tmp/epidemic',w=8,h=4)
# plot.network(E$P$G,list(fill='health.tf'),list()); fig.save('.tmp/network',w=8,h=6)

# build + run many
N.s = 7
P.s = def.params.s(N.s)
E.s = epi.run.s(P.s)
out.long.s = epi.output.melt.s(E.s)

# plot prevalence
g = plot.epidemic(out.long.s,select=list(var='N',health=c('S','E','I','H','R','V1'))) +
  labs(y='Count',color='State',fill='State') +
  facet_wrap('~health',scales='free_y')
  fig.save('.tmp/health',w=10,h=6)
g = plot.epidemic(out.long.s,select=list(var='inc',health='all')) +
  labs(y='Incidence',color='State',fill='State')
  fig.save('.tmp/inc',w=8,h=4)

q()
# plot networks
par.lapply(seq(N.s),function(s){
  plot.network(E.s[[s]]$P$G,list(fill='health.tf',color='inf.src',shape='inf.src'))
  fig.save('.tmp/net-',s,w=7,h=7)
})
