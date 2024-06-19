source('utils/all.r')
source('graph/graph.r')
source('model/meta.r')
source('model/params.r')
source('model/epidemic.r')

.debug = FALSE
P.args = list(seed=0,N.I0=10,N.V0=c(0,0))
mbm = microbenchmark::microbenchmark(
  {R = epi.run(kw.call(def.params,P.args,N=100))},
  {R = epi.run(kw.call(def.params,P.args,N=1000))},
  {R = epi.run(kw.call(def.params,P.args,N=10000))},
  {R = epi.run(kw.call(def.params,P.args,N=100000))},
times=10)
tmp.out(mbm,'epi.run',path=c('.tmp','mbm'))
