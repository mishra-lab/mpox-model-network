source('utils/all.r')
source('model/params.r')
source('model/network.r')
source('model/epidemic.r')
source('model/plot.r')
library('ggplot2')

t = epi.t(tf=180)

# DEBUG: run one
# P = def.params(seed=0)
# R = epi.results(P,t,epi.run(P,t))
# fig.save(quoted=quote(draw.network(R$P$G,n.color='health')),'Rplots',w=7,h=7)

# build + run many
N.s = 8
P.s = def.params.s(N.s)
R.s = epi.run.s(P.s,t)
out.long.s = epi.output.melt.s(R.s,P.s)

# plot prevalence by city
g = plot.epidemic(out.long.s,y='n.city',select=list(city=c('A','B','C')),facet='~city') +
  labs(y='Prevalence',color='State',fill='State')
fig.save('.tmp/prev.city',w=8,h=4)

# plot networks
for (s in seq(N.s)){
  fig.save(paste0('.tmp/net-',s),quote(draw.network(R.s[[s]]$P$G,n.color='health',shape='src')))
}