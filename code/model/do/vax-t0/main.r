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
  R.s = epi.run.s(P.s,t.vec)
  out.long.s = epi.output.melt.s(R.s)
}

vt0.grid.cia = function(N.s,vax.cov.v,vax.eff.v){
  # base case: no vaccine
  out.long.v0 = row.select(vt0.run(N.s,0,0),var='inc',health='all')
  cum.inf.v0 = cumsum.len(out.long.v0$value,n=N.s)
  # grid of cov & eff
  out.long =
    do.call(rbind,lapply(vax.eff.v,function(vax.eff){
      do.call(rbind,lapply(vax.cov.v,function(vax.cov){
        if (vax.eff==0 || vax.cov==0){ # don't run dummy cases
          out.long.v = out.long.v0
        } else {
          out.long.v = row.select(vt0.run(N.s,vax.cov,vax.eff),var='inc',health='all')
        }
        out.long.v$var = 'cia'
        out.long.v$value = 100 * (cum.inf.v0 - cumsum.len(out.long.v$value,,n=N.s)) / cum.inf.v0
        out.long.v$vax.cov = vax.cov
        out.long.v$vax.eff = vax.eff
        return(out.long.v)
      }))
  }))
}

vt0.plot = function(N.s=1,vax.cov.v=c(.2,.4,.6,.8),vax.eff.v=c(.2,.4,.6,.8)){
  # run the grid
  out.long = vt0.grid.cia(N.s=N.s,vax.cov.v=vax.cov.v,vax.eff.v=vax.eff.v)
  out.long$vax.cov = factor(out.long$vax.cov,vax.cov.v,paste0(100*vax.cov.v,'% Coverage'))
  out.long$vax.eff = factor(out.long$vax.eff,vax.eff.v,paste0(100*vax.eff.v,'% Effect.'))
  # plot
  g = plot.epidemic(out.long,select=list(var='cia'),color='vax.eff',conf.int=.5) +  
    facet_wrap('vax.cov') +
    scale_color_viridis(discrete=TRUE) + scale_fill_viridis(discrete=TRUE) +
    labs(color='',fill='',y='Cumulative Infections Averted (%)')
  fig.save('vt0-cia-plot',w=8,h=6)
}

vt0.surface = function(N.s=10,N.grid=7){
  # run the grid
  vax.cov.v = seq(0,1,l=N.grid)
  vax.eff.v = seq(0,1,l=N.grid)
  out.long = vt0.grid.cia(N.s=N.s,vax.cov.v=vax.cov.v,vax.eff.v=vax.eff.v)
  out.long.agg = aggregate(formula('value ~ vax.cov + vax.eff'),row.select(out.long,t=tf-1),median)
  # plot
  g = ggplot(out.long.agg,aes_string(x='100*vax.cov',y='100*vax.eff',z='value')) +
    geom_contour_filled(alpha=.8) +
    theme_light() +
    labs(x='Coverage (%)',y='Effect. (%)',fill='CIA (%)')
  fig.save('vt0-cia-surface',w=6,h=5)
}

vt0.plot()
vt0.surface()



