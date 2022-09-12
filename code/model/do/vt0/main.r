source('utils/all.r')
source('graph/graph.r')
source('model/meta.r')
source('model/params.r')
source('model/epidemic.r')
source('model/plot.r')

.debug = FALSE
vt0 = list()
vt0$tf = 180
# base
vt0$base$N.s     = 100
vt0$base$N       = 10000
vt0$base$vax.cov = .50
vt0$base$vax.eff = .50
# test
vt0$test$N.s     = 7
vt0$test$N       = 1000
vt0$test$vax.cov = c(.33,.67)
vt0$test$vax.eff = c(.33,.67)
# objective 1
vt0$obj1$N.s     = 10
vt0$obj1$N       = 10000
vt0$obj1$vax.cov = c(.2,.4,.6,.8)
vt0$obj1$vax.eff = c(.2,.4,.6,.8)
# objective 2a
vt0$obj2a$N.s       = 10
vt0$obj2a$N         = c(1000,3000,10000,30000) # + 1e5 ?
vt0$obj2a$par.scale = c(.5,.75,1,1.25,1.5)
vt0$obj2a$vax.cov   = .50
vt0$obj2a$vax.eff   = .50
# objective 2b
vt0$obj2b$N.s          = 10
vt0$obj2b$N            = 10000
vt0$obj2b$dur.IH.scale = c(.5,.75,1,1.25,1.5)
vt0$obj2b$p.asymp      = c(.05,.10,.15,.20,.25)
vt0$obj2b$vax.cov      = .50
vt0$obj2b$vax.eff      = .50

vt0.fname.fig   = function(slug){ root.path('out','fig','vt0',slug,create=TRUE) }
vt0.fname.rdata = function(slug){ root.path('data','.rdata','vt0',paste0(slug,'.rdata'),create=TRUE) }

vt0.run = function(args){
  # clean params
  print(paste(names(args),' = ',args,collapse=', '))
  P.s = kw.call(def.params.s,
    list.update(args,
      seeds        = 1:args$N.s,
      N.V0         = c(round(args$N*args$vax.cov),0),
      vax.eff.dose = c(args$vax.eff,0)))
  # run the model
  t.vec = epi.t(tf=vt0$tf)
  E.s = epi.run.s(P.s,t.vec)
  # clean / compute outputs
  out.long.s.all = epi.output.melt.s(E.s)
  out.long.s = cbind(rbind(
      row.select(out.long.s.all,var='inc',health='all'),
      row.select(out.long.s.all,var='inc',health='S'),
      row.select(out.long.s.all,var='inc',health='V1'),
      row.select(out.long.s.all,var='prev',health='E'),
      row.select(out.long.s.all,var='prev',health='I'),
      row.select(out.long.s.all,var='prev',health='H'),
      row.select(out.long.s.all,var='prev',health='R'),
      data.frame(t=NA,value=sapply(E.s,epi.tex),var='tex',health='all',seed=1:args$N.s)),
    args)
}

vt0.run.grid = function(args,out.i0=FALSE){
  # do vt0.run for a range of vax & other params, with vt0[[base]] as default params,
  # comparing vax > 0 vs vax=0 within each set of other params, and for multiple seeds
  args.vax   = expand.grid(list(vax.cov=args$vax.cov,vax.eff=args$vax.eff))
  args.other = expand.grid(list.update(args,vax.cov=NULL,vax.eff=NULL))
  out.long = do.call(rbind,lapply(1:nrow(args.other),function(i){
    out.long.i0 = vt0.run(as.list(cbind(args.other[i,],vax.cov=0,vax.eff=0))) # run no vax case
    b.inc = out.long.i0$var=='inc' # pre-compute inc selector
    cum.inf.fun = function(out.long.v){ chunk.fun(out.long.v[b.inc,]$value,l=vt0$tf,fun=cumsum) }
    cum.inf.i0 = cum.inf.fun(out.long.i0) # cumulative infections
    out.long.i0 = rbind(out.long.i0,list.update(out.long.i0[b.inc,],var='ci',value=cum.inf.i0))
    if (!out.i0){ out.long.i0 = NULL } # include out.long.i0 in output?
    out.long.i = rbind(out.long.i0,do.call(rbind,lapply(1:nrow(args.vax),function(v){ # vax grid
      out.long.iv = vt0.run(as.list(cbind(args.other[i,],args.vax[v,]))) # run this vax case
      # cumulative infections averted
      cum.inf.iv  = cum.inf.fun(out.long.iv)
      out.long.iv = rbind(out.long.iv,
        list.update(out.long.iv[b.inc,],var='ci', value=cum.inf.iv),
        list.update(out.long.iv[b.inc,],var='cia',value=100*(cum.inf.i0-cum.inf.iv)/cum.inf.i0)
      )
    })))
  }))
}

vt0.plot.tex = function(out.long,facet=NULL,...){
  # see https://stackoverflow.com/questions/73480520
  N.s = len(unique(out.long$seed))
  # TODO: compute group automatically?
  noise = function(x){ x + runif(len(x),-1e-9,+1e-9) } # avoid unique
  cdf.adj = function(y,...){ unlist(lapply(split(y,list(...)),function(yi){
    yi * (len(yi)-2) / N.s # only non-NA in yi, so we use len(yi) vs N.s to adjust; -2 from pad
  })) }
  g = ggplot(row.select(out.long,var='tex'),aes_string(x='noise(value)',...)) +
    stat_ecdf(aes(y=after_stat(cdf.adj(y,group,PANEL)))) +
    facet_grid(facet) +
    lims(x=c(0,vt0$tf),y=c(0,1)) +
    labs(y=.labs$tex,x='Time to Extinction (days)')
  g = plot.clean(g)
}

vt0.clean.out.long = function(out.long,args){
  b = out.long$var=='inc'
  out.long[b,]$value = out.long[b,]$value / out.long[b,]$N * 365
  b = out.long$var=='prev'
  out.long[b,]$value = out.long[b,]$value * 100
  f = function(name,scale=1){
    if (is.null(out.long[[name]])){ return(NULL) }
    factor(out.long[[name]],args[[name]],paste0(scale*args[[name]])) }
  out.long$vax.cov      = f('vax.cov',100)
  out.long$vax.eff      = f('vax.eff',100)
  out.long$N            = f('N',1)
  out.long$par.scale    = f('par.scale',1)
  out.long$p.asymp      = f('p.asymp',100)
  out.long$dur.IH.scale = f('dur.IH.scale',1)
  return(out.long)
}

.labs = list(
  'scen'         = 'Scenario',
  'vax.cov'      = 'Coverage (%)',
  'vax.eff'      = 'Effectiveness (%)',
  'N'            = 'Population Size',
  'par.scale'    = 'Relative\nPartner Numbers',
  'p.asymp'      = 'Proportion Asymptomatic (%)',
  'dur.IH.scale' = 'Relative\nTime to Isolation',
  'ci'           = 'Cumulative Infections',
  'cia'          = 'Cumulative Infections Averted (%)',
  'tex'          = 'Cumulative Probability of Extinction (%)',
  'inc'          = 'Incidence (per person-year)',
  'prev'         = 'Infection Prevalence (%)')

vt0.example.cia = function(){
  vars = c('Absolute Incidence'='inc','Cumulative Infections'='ci','Infections Averted (%)'='cia')
  out.long = vt0.run.grid(list.update(vt0$base),out.i0=TRUE) # run base & vax scenarios
  out.long = row.select(out.long,var=vars,health='all') # select outputs
  out.long$varf = factor(out.long$var,levels=vars,labels=names(vars)) # pretty factor names
  out.long$scen = factor(out.long$vax.cov,labels=c('No Vaccine','50% Coverage')) # pretty color names
  g = plot.epidemic(out.long,select=list(var=vars),color='scen',facet='varf',scales='free_y') +
    scale_color_manual(values=as.vector(M$health$color[c('I','V1')])) +
    scale_fill_manual(values=as.vector(M$health$color[c('I','V1')])) +
    labs(y='Cases',color=.labs$scen,fill=.labs$scen) +
    theme(legend.position='top')
  fig.save(vt0.fname.fig('example-cia'),w=5,h=7)
}

vt0.obj = function(case,.run=FALSE){
  args = vt0[[case]]
  slug = paste0(case,'-s',args$N.s)
  if (.run){
    out.long = vt0.run.grid(args)
    save(file=vt0.fname.rdata(slug),out.long)
  } else {
    load(file=vt0.fname.rdata(slug))
  }
  out.long = vt0.clean.out.long(out.long,args)
  plot.save = function(g,which){
    glab = .labs[[plot.args$color]]
    g = g + labs(color=glab,fill=glab) +
      scale_color_viridis(drop=FALSE,discrete=is.factor(out.long[[plot.args$color]])) +
      scale_fill_viridis(drop=FALSE,discrete=is.factor(out.long[[plot.args$color]])) +
      theme(legend.position='top')
    fig.save(vt0.fname.fig(paste0(slug,'-',which)),w=5,h=8)
  }
  plot.args = c(list(out.long=out.long),list(
      'test'  = list(facet='vax.eff',color='vax.cov'),
      'obj1'  = list(facet='vax.eff',color='vax.cov'),
      'obj2a' = list(facet='N',      color='par.scale'),
      'obj2b' = list(facet='p.asymp',color='dur.IH.scale')
    )[[case]])
  g = kw.call(vt0.plot.tex,plot.args,group=plot.args$color); plot.save(g,'tex')
  g = kw.call(plot.epidemic,plot.args,select=list(var='cia',health='all')) + labs(y=.labs$cia); plot.save(g,'cia')
  g = kw.call(plot.epidemic,plot.args,select=list(var='inc',health='all')) + labs(y=.labs$inc); plot.save(g,'inc')
  g = kw.call(plot.epidemic,plot.args,select=list(var='prev',health='I')) + labs(y=.labs$prev); plot.save(g,'prev')
}

if (sys.nframe() == 0){
  vt0.example.cia()
  # vt0.obj('obj1')
  # vt0.obj('obj2a')
  # vt0.obj('obj2b')
}
