source('utils/all.r')
source('graph/graph.r')
source('model/meta.r')
source('model/params.r')
source('model/epidemic.r')
source('model/plot.r')

.debug = FALSE

t = epi.t(tf=180)

# DEBUG: run one
# P = def.params(seed=0)
# R = epi.run(P,t)
# g = plot.epidemic(epi.output.melt(R$out,P)); fig.save('.tmp/epidemic',w=8,h=4)
# plot.network(R$P$G,list(fill='health'),list()); fig.save('.tmp/network',w=8,h=6)

# build + run many
N.s = 8
P.s = def.params.s(N.s)
R.s = epi.run.s(P.s,t)
out.long.s = epi.output.melt.s(R.s,P.s)

# plot prevalence
g = plot.epidemic(out.long.s,select=list(var='N',health=c('S','E','I','R'))) +
  labs(y='Count',color='State',fill='State') +
  facet_wrap('~health',scales='free_y')
  print(g); q()
  fig.save('.tmp/health',w=10,h=6,ext='png')
g = plot.epidemic(out.long.s,select=list(var='inc',health='all')) +
  labs(y='Incidence',color='State',fill='State')
  fig.save('.tmp/inc',w=8,h=4,ext='png')

q()
# plot networks
par.lapply(seq(N.s),function(s){
  plot.network(R.s[[s]]$P$G,list(fill='health',color='inf.src',shape='inf.src'))
  fig.save('.tmp/net-',s,w=7,h=7)
})
