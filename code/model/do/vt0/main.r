source('utils/all.r')
source('graph/graph.r')
source('model/meta.r')
source('model/params.r')
source('model/epidemic.r')
source('model/plot.r')

# config - TODO: surface stuff here too?
.debug = FALSE
vt0 = list()
vt0$tf = 360
vt0$N.s = 10
vt0$N = 1000
vt0$N.v = c(1000,3000,10000,30000)
vt0$vax.cov.v = c(.2,.4,.6,.8)
vt0$vax.eff.v = c(.2,.4,.6,.8)

.fname = function(slug,ext,...){
  f = paste0('vt0-',slug,'-t',vt0$tf,'-s',vt0$N.s,'-v',len(vt0$N.v),ext)
  root.path(...,'vt0',f,create=TRUE)
}

vt0.fname.fig = function(slug){ .fname(slug,'','out','fig') }
vt0.fname.rdata = function(slug){ .fname(slug,'.rdata','data','.rdata') }

.clean.out.long = function(out.long,N.v,vax.cov.v,vax.eff.v){
  b.inc = out.long$var=='inc'
  out.long[b.inc,]$value = out.long[b.inc,]$value / out.long[b.inc,]$N
  out.long$N       = factor(out.long$N,N.v,paste0('N = ',N.v))
  out.long$vax.cov = factor(out.long$vax.cov,vax.cov.v,paste0(100*vax.cov.v,'% Coverage'))
  out.long$vax.eff = factor(out.long$vax.eff,vax.eff.v,paste0(100*vax.eff.v,'% Effect.'))
  return(out.long)
}

.var.label = list(
  'cia'  = 'Cumulative Infections Averted (%)',
  'inc'  = 'Incidence',
  'prev' = 'Prevalence'
)

vt0.run = function(seeds,N,vax.cov,vax.eff){
  print(sprintf('N = %d, cov = %.2f, eff = %.2f',N,vax.cov,vax.eff))
  t.vec = epi.t(tf=vt0$tf)
  P.s = def.params.s(seeds,N=N,N.V0=c(round(N*vax.cov),0),vax.eff.dose=c(vax.eff,0))
  E.s = epi.run.s(P.s,t.vec)
  out.long.s.all = epi.output.melt.s(E.s)
  out.long.s = cbind(rbind(
      row.select(out.long.s.all,var='inc',health='all'),
      row.select(out.long.s.all,var='inc',health='S'),
      row.select(out.long.s.all,var='inc',health='V1'),
      row.select(out.long.s.all,var='prev',health='E'),
      row.select(out.long.s.all,var='prev',health='I'),
      row.select(out.long.s.all,var='prev',health='H'),
      row.select(out.long.s.all,var='prev',health='R'),
      data.frame(t=NA,value=sapply(E.s,epi.tex),var='tex',health='all',seed=seeds)), # extinction
    N=N,vax.cov=vax.cov,vax.eff=vax.eff)
}

vt0.run.grid = function(N.s,N.v,vax.cov.v,vax.eff.v){
  # base case: no vaccine
  out.long.v0 = lapply(N.v,vt0.run,seeds=1:N.s,vax.cov=0,vax.eff=0)
  b.inc = row.select(out.long.v0[[1]],var='inc',.return='b') # pre-compute inc selector
  b.tex = row.select(out.long.v0[[1]],var='tex',.return='b') # pre-compute tex selector
  cum.inf.fun = function(out.long.v){ chunk.fun(out.long.v[b.inc,]$value,n=N.s,fun=cumsum) }
  cum.inf.v0 = lapply(out.long.v0,cum.inf.fun) # cumulative infections
  # grid of cov & eff
  out.long =
  do.call(rbind,lapply(1:len(N.v),function(i){
    do.call(rbind,lapply(vax.eff.v,function(vax.eff){
      do.call(rbind,lapply(vax.cov.v,function(vax.cov){
        if (vax.eff==0 || vax.cov==0){ # don't run dummy cases
          out.long.v = out.long.v0[[i]]
        } else {
          out.long.v = vt0.run(N.s,N.v[i],vax.cov,vax.eff)
        }
        # cumulative infections averted
        out.long.v.cia = out.long.v[b.inc,]
        out.long.v.cia$var = 'cia'
        out.long.v.cia$value = 100 * (cum.inf.v0[[i]] - cum.inf.fun(out.long.v)) / cum.inf.v0[[i]]
        return(rbind(out.long.v,out.long.v.cia))
      }))
    }))
  }))
}

vt0.plot.var = function(out.long,var='cia',health='all',conf.int=.9){
  # plot
  select = list(var=var,health=health)
  g = plot.epidemic(out.long,select=select,color='vax.eff',facet='vax.cov ~ N',conf.int=conf.int) +
    scale_color_viridis(discrete=TRUE) + scale_fill_viridis(discrete=TRUE) +
    labs(color='',fill='',y=.var.label[[var]])
  fig.save(vt0.fname.fig(paste0('plot-',var)),w=10,h=8)
}

vt0.plot.tex = function(out.long){
  # stupid hack as ecdf does not support NA
  noise = function(x){ x + runif(len(x),-1e-9,+1e-9) } # avoid unique
  cdf.adj = function(y,...){ unlist(lapply(split(y,list(...)),function(yi){
    yi * (len(yi)-2) / vt0$N.s # only non-NA in yi, so we use len(yi) vs N.s to adjust; -2 from pad
  })) }
  g = ggplot(row.select(out.long,var='tex'),aes(x=noise(value),color=vax.eff,fill=vax.eff)) +
    stat_ecdf(aes(y=after_stat(cdf.adj(y,group,PANEL)))) +
    facet_grid('vax.cov ~ N') +
    theme_light() +
    lims(x=c(0,vt0$tf)) +
    scale_color_viridis(discrete=TRUE) + scale_fill_viridis(discrete=TRUE) +
    labs(color='',fill='',y='Cumulative Probability of Extinction',x='Time to Extinction (days)')
  fig.save(vt0.fname.fig('plot-tex'),w=10,h=8)
}

vt0.surface.cia = function(N.s=10,N.grid=7){
  # run the grid
  vax.cov.v = seq(0,1,l=N.grid)
  vax.eff.v = seq(0,1,l=N.grid)
  out.long = vt0.run.grid(N.s,vt0$N,vax.cov.v,vax.eff.v)
  out.long = row.select(out.long,var='cia',t=vt0$tf-1)
  out.long.agg = aggregate(formula('value ~ vax.cov + vax.eff'),out.long,median)
  # plot
  g = ggplot(out.long.agg,aes_string(x='100*vax.cov',y='100*vax.eff',z='value')) +
    geom_contour_filled(alpha=.8) +
    theme_light() +
    labs(x='Coverage (%)',y='Effect. (%)',fill='CIA (%)')
  fig.save(vt0.fname.fig('surface-cia'),w=10,h=8)
}

if (sys.nframe() == 0){
  # out.long = do.call(vt0.run.grid,vt0[c('N.s','N.v','vax.cov.v','vax.eff.v')])
  # save(file=vt0.fname.rdata('out-long'),out.long); q()
  load(file=vt0.fname.rdata('out-long'))
  out.long = .clean.out.long(out.long,vt0$N.v,vt0$vax.cov.v,vt0$vax.eff.v)
  par.funs(out.long=out.long,list( # generate plots in parallel
    list(vt0.plot.var,var='cia'),
    list(vt0.plot.var,var='inc',health='all'),
    list(vt0.plot.var,var='prev',health='I'),
    list(vt0.plot.tex)
  ))
  # vt0.surface.cia()
}
