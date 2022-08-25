# load
source('model/do/vax-t0/main.r')
# config
.n.cores = 80
vt0$N.s = 1000
# run + plot
out.long = do.call(vt0.run.grid,vt0[c('N.s','N.v','vax.cov.v','vax.eff.v')])
save(file=vt0.fname.rdata('out-long'),out.long)
out.long = .clean.out.long(out.long,vt0$N.v,vt0$vax.cov.v,vt0$vax.eff.v)
vt0.plot.cia(out.long)
vt0.plot.tex(out.long)