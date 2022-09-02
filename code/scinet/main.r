# load
source('model/do/vt0/main.r')
# config
.n.cores = 80
vt0$N.s = 100
# run + plot
out.long = do.call(vt0.run.grid,vt0[c('N.s','N.v','vax.cov.v','vax.eff.v')])
save(file=vt0.fname.rdata('out-long'),out.long)
load(file=vt0.fname.rdata('out-long'))
out.long = .clean.out.long(out.long,vt0$N.v,vt0$vax.cov.v,vt0$vax.eff.v)
par.funs(out.long=out.long,list( # generate plots in parallel
  list(vt0.plot.var,var='cia'),
  list(vt0.plot.var,var='inc',health='all'),
  list(vt0.plot.var,var='prev',health='I'),
  list(vt0.plot.tex)
))
