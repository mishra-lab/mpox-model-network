source('utils/all.r')
source('graph/graph.r')
source('model/meta.r')
source('model/params.r')
source('model/epidemic.r')
source('model/plot.r')

.debug = FALSE
tf = 180
N  = 1000

vt0.run = function(seeds,vax.cov,vax.eff){
  print(sprintf('cov = %.2f, eff = %.2f',vax.cov,vax.eff))
  t.vec = epi.t(tf=tf)
  P.s = def.params.s(seeds,N=N,N.V0=c(round(N*vax.cov),0),vax.eff.dose=c(vax.eff,0))
  E.s = epi.run.s(P.s,t.vec)
  out.long.s = cbind(rbind(
      row.select(epi.output.melt.s(E.s),var='inc',health='all'), # incidence
      data.frame(t=NA,value=sapply(E.s,epi.tex),var='tex',health='all',seed=seeds)), # extinction
    vax.cov=vax.cov,vax.eff=vax.eff)
}

vt0.run.grid = function(N.s,vax.cov.v,vax.eff.v){
  # base case: no vaccine
  out.long.v0 = vt0.run(1:N.s,0,0)
  b.inc = row.select(out.long.v0,var='inc',.return='b') # pre-compute inc selector
  b.tex = row.select(out.long.v0,var='tex',.return='b') # pre-compute tex selector
  cum.inf.fun = function(out.long.v){ chunk.fun(out.long.v[b.inc,]$value,n=N.s,fun=cumsum) }
  cum.inf.v0 = cum.inf.fun(out.long.v0) # cumulative infections
  # grid of cov & eff
  out.long =
    do.call(rbind,lapply(vax.eff.v,function(vax.eff){
      do.call(rbind,lapply(vax.cov.v,function(vax.cov){
        if (vax.eff==0 || vax.cov==0){ # don't run dummy cases
          out.long.v = out.long.v0
        } else {
          out.long.v = vt0.run(N.s,vax.cov,vax.eff)
        }
        # cumulative infections averted
        out.long.v[b.inc,]$var = 'cia'
        out.long.v[b.inc,]$value = 100 * (cum.inf.v0 - cum.inf.fun(out.long.v)) / cum.inf.v0
        return(out.long.v)
      }))
  }))
}

vt0.plot.cia = function(N.s=1,vax.cov.v=c(.2,.4,.6,.8),vax.eff.v=c(.2,.4,.6,.8),conf.int=.5){
  # run the grid
  out.long = row.select(vt0.run.grid(N.s=N.s,vax.cov.v=vax.cov.v,vax.eff.v=vax.eff.v),var='cia')
  out.long$vax.cov = factor(out.long$vax.cov,vax.cov.v,paste0(100*vax.cov.v,'% Coverage'))
  out.long$vax.eff = factor(out.long$vax.eff,vax.eff.v,paste0(100*vax.eff.v,'% Effect.'))
  # plot
  g = plot.epidemic(out.long,select=list(var='cia'),color='vax.eff',conf.int=conf.int) +
    facet_wrap('vax.cov') +
    scale_color_viridis(discrete=TRUE) + scale_fill_viridis(discrete=TRUE) +
    labs(color='',fill='',y='Cumulative Infections Averted (%)')
  fig.save('vt0-plot-cia',w=8,h=6)
}

vt0.plot.tex = function(N.s=1,vax.cov.v=c(.2,.4,.6,.8),vax.eff.v=c(.2,.4,.6,.8)){
  # run the grid
  out.long = row.select(vt0.run.grid(N.s=N.s,vax.cov.v=vax.cov.v,vax.eff.v=vax.eff.v),var='tex')
  out.long$vax.cov = factor(out.long$vax.cov,vax.cov.v,paste0(100*vax.cov.v,'% Coverage'))
  out.long$vax.eff = factor(out.long$vax.eff,vax.eff.v,paste0(100*vax.eff.v,'% Effect.'))
  # plot
  g = ggplot(out.long,aes_string(x='value',y='vax.eff',color='vax.eff',fill='vax.eff')) +
    geom_density_ridges(alpha=.5) +
    facet_wrap('vax.cov') +
    theme_light() +
    scale_color_viridis(discrete=TRUE) + scale_fill_viridis(discrete=TRUE) +
    labs(color='',fill='',y='',x='Time to Extinction (days)')
  fig.save('vt0-plot-tex',w=8,h=6)
}

vt0.surface.cia = function(N.s=10,N.grid=7){
  # run the grid
  vax.cov.v = seq(0,1,l=N.grid)
  vax.eff.v = seq(0,1,l=N.grid)
  out.long = row.select(vt0.run.grid(N.s=N.s,vax.cov.v=vax.cov.v,vax.eff.v=vax.eff.v),var='cia',t=tf-1)
  out.long.agg = aggregate(formula('value ~ vax.cov + vax.eff'),out.long,median)
  # plot
  g = ggplot(out.long.agg,aes_string(x='100*vax.cov',y='100*vax.eff',z='value')) +
    geom_contour_filled(alpha=.8) +
    theme_light() +
    labs(x='Coverage (%)',y='Effect. (%)',fill='CIA (%)')
  fig.save('vt0-surface-cia',w=6,h=5)
}

vt0.plot.cia(N.s=10,conf.int=.9)
vt0.plot.tex(N.s=10)
vt0.surface.cia()



