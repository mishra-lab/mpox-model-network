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
vt0$N.s = 100
vt0$N = 1000
vt0$N.v = c(1000,3000,10000,30000)
vt0$vax.cov.v = c(.2,.4,.6,.8)
vt0$vax.eff.v = c(.2,.4,.6,.8)

.clean.out.long = function(out.long,N.v,vax.cov.v,vax.eff.v){
  out.long$N       = factor(out.long$N,N.v,paste0(N.v))
  out.long$vax.cov = factor(out.long$vax.cov,vax.cov.v,paste0(100*vax.cov.v,'% Coverage'))
  out.long$vax.eff = factor(out.long$vax.eff,vax.eff.v,paste0(100*vax.eff.v,'% Effect.'))
  return(out.long)
}

vt0.run = function(seeds,N,vax.cov,vax.eff){
  print(sprintf('N = %d, cov = %.2f, eff = %.2f',N,vax.cov,vax.eff))
  t.vec = epi.t(tf=vt0$tf)
  P.s = def.params.s(seeds,N=N,N.V0=c(round(N*vax.cov),0),vax.eff.dose=c(vax.eff,0))
  E.s = epi.run.s(P.s,t.vec)
  out.long.s = cbind(rbind(
      row.select(epi.output.melt.s(E.s),var='inc',health='all'), # incidence
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
        out.long.v[b.inc,]$var = 'cia'
        out.long.v[b.inc,]$value = 100 * (cum.inf.v0[[i]] - cum.inf.fun(out.long.v)) / cum.inf.v0[[i]]
        return(out.long.v)
      }))
    }))
  }))
}

vt0.plot.cia = function(out.long,conf.int=.9){
  # plot
  g = plot.epidemic(out.long,select=list(var='cia'),color='vax.eff',conf.int=conf.int) +
    facet_grid('vax.cov ~ N') +
    scale_color_viridis(discrete=TRUE) + scale_fill_viridis(discrete=TRUE) +
    labs(color='',fill='',y='Cumulative Infections Averted (%)')
  fig.save('vt0-plot-cia',w=8,h=6)
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
  fig.save('vt0-plot-tex',w=10,h=8)
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
  fig.save('vt0-surface-cia',w=10,h=8)
}
if (sys.nframe() == 0){
  out.long = do.call(vt0.run.grid,vt0[c('N.s','N.v','vax.cov.v','vax.eff.v')])
  out.long = .clean.out.long(out.long,vt0$N.v,vt0$vax.cov.v,vt0$vax.eff.v)
  vt0.plot.cia(out.long)
  vt0.plot.tex(out.long)
  vt0.surface.cia()
}
