source('utils/all.r')
source('graph/graph.r')
source('model/meta.r')
source('model/params.r')
source('model/epidemic.r')
source('model/plot.r')

.debug = FALSE
t = epi.t(tf=180)
P.args = list(seed=0,N.I0=10,N.V0=c(0,0))
mbm = microbenchmark::microbenchmark(
  {R = epi.run(kw.call(def.params,P.args,N=100),t)},
  {R = epi.run(kw.call(def.params,P.args,N=1000),t)},
  {R = epi.run(kw.call(def.params,P.args,N=10000),t)},
  {R = epi.run(kw.call(def.params,P.args,N=100000),t)},
times=10)
tmp.out(mbm,'epi.run',path=c('.tmp','mbm'))