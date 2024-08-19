# a collection of plots & tools to hand fit the epidemic

source('utils/all.r')
source('graph/graph.r')
source('model/meta.r')
source('model/params.r')
source('model/epidemic.r')
source('model/plot.r')

fname = function(slug){
  # e.g. .../out/fig/.tmp/fit/2020-01-01/{slug}
  root.path('out','fig','.tmp','fit',Sys.Date(),slug,create=TRUE)
}

# -----------------------------------------------------------------------------
# num ptrs vs data

x.max.dt = function(dt){ list('7'=10,'30'=20,'180'=50)[[as.character(dt)]] }

fit.ptr.engage = function(N.s=21){
  source('data/engage.r')
  X.dat = row.select(main.p6m(main.stat()),city='avg')
  fit = cbind(#    init,   min,  max
    w.sdlog    = c( 1.3,   0.4  , 9.0),
    r.ptr.casu = c( 0.020, 0.001, 0.3),
    r.ptr.once = c( 0.020, 0.001, 0.3),
    w.pwr.excl = c(+0.200,-1    ,+1  ),
    w.pwr.open = c(-0.300,-1    ,+1  ))
  fit.fix = list(
    d.main.scale = 1059,
    d.casu.scale = 299,
    p.excl =  0.154,
    p.open =  0.309,
    p.adj  =  0)
  err.fun = function(fit){ cat(fit)
    X.mod = fit.ptr.mod(N.s=N.s,fit=c(as.list(fit),fit.fix),status='main.now')
    e = sum(abs(X.dat$cp - X.mod$cp)^2)/N.s
    cat(' > e = ',e,'\n'); return(e) }
  # opt = optim(unlist(fit[1,]),err.fun,method='L-BFGS-B',
  #   lower=unlist(fit[2,]),upper=unlist(fit[3,]))
  opt = list(par=unlist(fit[1,])); opt$value = err.fun(opt$par) # DEBUG
  print(opt)
  X.mod = fit.ptr.mod(N.s=100,fit=c(as.list(opt$par),fit.fix),status='main.now')
  fit.ptr.plot(X.dat,X.mod); fig.save(fname('ptr.180.engage'),w=5,h=5)
}

fit.ptr.kenya = function(N.s=21){
  source('data/kenya.r')
  X.dat = subset(main.ptr(),x < 6 | county=='all')
  fit = cbind(
    w.sdlog      = c( 0.4  , 0.1  , 10 ),
    r.ptr.casu   = c( 0.005, 0.001, 0.1),
    r.ptr.once   = c( 0.030, 0.001, 0.1),
    w.pwr.excl   = c( 0.000,-3    ,+3  ),
    w.pwr.open   = c( 0.000,-3    ,+3  ),
    d.main.scale = c(  720 , 30   ,3600),
    d.casu.scale = c(   60 , 30   ,3600),
    p.excl       = c( 0.100, 0.001,0.35),
    p.open       = c( 0.500, 0.001,0.75),
    p.adj        = c( 0.200, 0.001,0.50))
  X.dat.7  = subset(X.dat,county=='all' & dt== 7)
  X.dat.30 = subset(X.dat,county=='all' & dt==30)
  err.fun = function(fit){ cat(fit)
    e = sum(c(
      (X.dat.7$cp  - fit.ptr.mod(N.s=N.s,fit=as.list(fit),dt= 7,x=0:20)$cp)^2,
      (X.dat.30$cp - fit.ptr.mod(N.s=N.s,fit=as.list(fit),dt=30,x=0:20)$cp)^2))/N.s
    cat(' > e = ',e,'\n'); return(e) }
  # opt = optim(unlist(fit[1,]),err.fun,method='L-BFGS-B',
  #   lower=unlist(fit[2,]),upper=unlist(fit[3,]))
  opt = list(par=unlist(fit[1,])); opt$value = err.fun(opt$par) # DEBUG
  print(opt)
  for (dti in c(7,30)){
    x.max = x.max.dt(dti)
    X.mod = cbind(fit.ptr.mod(N.s=100,fit=as.list(opt$par),dt=dti,x=0:x.max),county='mod')
    X.dat.dt = subset(X.dat,dt==dti & x<=x.max)
    fit.ptr.plot(X.dat.dt,X.mod,dti,x.max=x.max,clr='county')
    fig.save(fname(sprintf('ptr.%d.kenya',dti)),w=5,h=5)
  }
}

fit.ptr.context = function(N.s=7){
  for (dti in c(30,180)){
    x.max = x.max.dt(dti)
    X = rbind(cbind(fit.ptr.mod(N.s,dt=dti,x=0:x.max,context='engage'),context='engage'),
              cbind(fit.ptr.mod(N.s,dt=dti,x=0:x.max,context='kenya'), context='kenya'))
    fit.ptr.plot(X[,],X,dti,x.max=x.max,clr='context')
    fig.save(fname(sprintf('ptr.%d.context',dti)),w=5,h=5)
  }
}

fit.ptr.mod = function(N.s,x=0:300,dt=180,...,N=1000,status='.'){
  # lookup or compute histogram of # partners observed in dt
  # for N.s param sets (pop size = N), possibly stratified by status
  t.obs = seq(0,180,max(dt,30))
  P.s = def.params.s(N.s,N=N,t.max=180,...,.par=FALSE)
  X.mod = rbind.lapply(P.s,function(P){ # each param set
    P$G$attr$i$. = 0 # dummy stat
    t0.e = P$G$attr$e$t0 # for speed
    tf.e = P$G$attr$e$tf # for speed
    i.stat = split(P$G$i,P$G$attr$i[[status]]) # split by status
    X.mod.s = rbind.lapply(names(i.stat),function(stat){ # each status
      i.s = i.stat[[stat]] # individuals with this status
      if (dt==P$t.max){ # lookup # ptrs -> 1 histogram
        nt = tabulate0(P$G$attr$i$n.ptr.tot[i.s],max(x))
        nt = matrix(nt,nrow=len(x)) # for apply below
      } else { # compute # ptrs -> len(t.obs) histograms
        e.i = lapply(i.s,function(i){ # partners per individual
          which(P$G$ii.e[,1]==i | P$G$ii.e[,2]==i) })
        p.adj = ifelse(is.null(P$fit$p.adj),0,P$fit$p.adj)
        a.i = runif(e.i) < p.adj # inds who report sex acts > num ptrs
        nt = sapply(t.obs,function(t){ # for each recall window
          tabulate0(mapply(function(e,a){ # histogram: reported ptrs
            if (a){ # num total sex acts (over-report)
              d.e = pmax(0, pmin(tf.e[e],t+dt) - pmax(t0.e[e],t))
              n = rbinom(1,sum(d.e),P$f.sex)
            } else { # num pts (ideal report)
              n = sum(t0.e[e] <= t+dt & tf.e[e] >= t)
            }
          },e.i,a.i),max(x))
        })
      }
      X.mod.ms = data.frame(seed=P$seed,dt=dt,status=stat,x=x,
        p=c(apply(nt,2,sum1)),cp=c(apply(nt,2,cum1)))
    })
  })
}

fit.ptr.plot = function(X.dat,X.mod,dt=180,x.max=50,clr='status'){
  ptr.facet = c(p='Proportion',cp='Cumulative Proportion')
  cols = intersect(names(X.dat),names(X.mod))
  X = rbind(cbind(X.dat[,cols],src='Data'),cbind(X.mod[,cols],src='Model'))
  X = melt(X,measure=names(ptr.facet),var='facet')
  X$facet = factor(X$facet,names(ptr.facet),ptr.facet)
  g = ggplot(X,aes(x=x,y=value,color=.data[[clr]],fill=.data[[clr]],linetype=src)) +
    facet_grid('facet',scales='free_y') +
    stat_summary(geom='ribbon',fun.min=min,fun.max=max,alpha=.3,color=NA) +
    stat_summary(geom='line',fun=median) +
    lims(x=c(0,x.max),y=c(0,NA)) +
    labs(y='Proportion',linetype='') +
    # scale_x_continuous(trans='log10') + # DEBUG
    xlab(sprintf('Reported number of partners in %d days',dt))
  g = add.meta.scales(g,list(color=clr,fill=clr))
  g = plot.clean(g)
}

# -----------------------------------------------------------------------------
# ptr durs vs data

fit.pdur.kenya = function(N.s=21,N.ptr=3,d.max=3650){
  source('data/kenya.r')
  x = 0:d.max
  P.s = def.params.s(N.s,N=1000,context='kenya')
  X.mod = rbind.lapply(P.s,function(P){
    t0.e = P$G$attr$e$t0 # speed
    tf.e = pmin(P$G$attr$e$tf,P$t.max) # speed + censoring
    once.e = P$G$attr$e$type == 'once' # speed
    dur.s = rbind.lapply(P$G$i,function(i){ # for each individual
      b.all  = P$G$ii.e[,1]==i | P$G$ii.e[,2]==i # ptrs of i
      e.all  = which(b.all) # all ptrs
      e.sub  = which(b.all & !once.e) # all ptrs except once
      dur.si = cbind( # durs of last N.ptr
        all = sort(sample(tf.e[e.all]-t0.e[e.all]),decr=TRUE)[1:N.ptr],
        sub = sort(sample(tf.e[e.sub]-t0.e[e.sub]),decr=TRUE)[1:N.ptr])
    })
    X.mod.s = data.frame(x=x,cp=c(ecdf(dur.s[,1])(x),ecdf(dur.s[,2])(x)),
      ptrs=rep(c('All','Non-Once'),each=len(x)))
  })
  X.dat = data.frame(x=x,ptrs=NA,cp=do.call(pweibull,c(list(q=x),args.pdur)))
  g = fit.pdur.plot(X.mod,X.dat,N.ptr)
  fig.save(fname('pdur.3r.kenya'),w=5,h=3)
}

fit.pdur.plot = function(X.mod,X.dat,N.ptr){
  X = rbind(cbind(X.mod,src='Sim'),cbind(X.dat,src='Data'))
  g = ggplot(X,aes(x=x,y=cp,color=ptrs,fill=ptrs,linetype=src)) +
    stat_summary(geom='ribbon',fun.min=min,fun.max=max,alpha=.2,color=NA) +
    stat_summary(geom='line',fun=median) +
    geom_vline(xintercept=30,linetype='dotted') + # 1 month ref
    scale_x_continuous(trans='log10') +
    labs(x=sprintf('Duration with last %d partners (days)',N.ptr),
         y='Cumulative Proportion',linetype='',
         color='Partners',fill='Partners') +
    ylim(c(0,1))
  g = plot.clean(g)
}

# -----------------------------------------------------------------------------
# other param plots

fit.durs.plot = function(){
  P = def.params()
  funs = list(
    dur.EI = 'Incubation Period',
    dur.IR = 'Infectious Period',
    dur.IH = 'Non-Isolated Period')
  clrs = unname(M$health$color[2:4])
  rfuns = setNames(P[paste0(names(funs),'.rfun')],funs)
  g = plot.durs(rfuns) +
    scale_color_manual(values=clrs) +
    scale_fill_manual(values=clrs)
    fig.save(fname('rfuns'),w=8,h=5)
}

# -----------------------------------------------------------------------------
# model outputs plots

fit.run.epi = function(seed=0,...){
  # run the model once
  P = def.params(seed=seed,...)
  E = epi.run(P)
}

fit.plot.epidemic = function(E){
  # main standard plots for one run
  # TODO: maybe add multiple model runs by default
  g = plot.epidemic(epi.output.melt(E),select=list(health=c('E','I','H'))) +
    labs(y='Individuals')
    fig.save(fname('epidemic'),w=8,h=4)
  g = plot.epidemic(epi.output.melt(E),select=list(var='inc')) +
    labs(y='Incidence')
    fig.save(fname('incidence'),w=8,h=4)
}

fit.plot.tree = function(E){
  # common stuff
  add.cmaps = function(g,label,...){ g = g +
    viridis::scale_color_viridis(...,discrete=TRUE) +
    viridis::scale_fill_viridis(...,discrete=TRUE) +
    labs(color=label,fill=label) }
  # transmission by generation
  E$tree$gen.f = factor(E$tree$gen)
  g = plot.tree(E$tree,fill='gen.f')
    g = add.tree.margin(add.cmaps(g,'Generation',option='inferno',begin=.2,end=.9))
    fig.save(fname('tree-gen'),g=g,w=8,h=6)
  # transmission by n children
  E$tree$n.chi.cut = int.cut(E$tree$n.chi.dir,c(0,1,2,6))
  g = plot.tree(E$tree,fill='n.chi.cut')
    g = add.tree.margin(add.cmaps(g,'Transmission',option='viridis'))
    fig.save(fname('tree-n-child'),g=g,w=8,h=6)
}

# -----------------------------------------------------------------------------
# dynamic network gif (expensive)

fit.network.gif = function(seed=0,N=100,tf=100,fps=5){
  # make a gif from plot.network, showing health states over time
  library('gganimate')
    # define params & run model
  P = def.params(seed=seed,N=N,N.I0=5)
  E = epi.run(P)
  t.vec = P$t.vec[1:tf]
  # hijack attr: write too many values that will be recycled by plot.network
  act.e = outer(t.vec,P$G$attr$e$t0,`>=`) & outer(t.vec,P$G$attr$e$tf,`<`)
  P$G$attr$i$t = rep(t.vec,each=P$G$N.i)
  P$G$attr$e$t = rep(t.vec,each=P$G$N.e)
  P$G$attr$e$active = factor(c(t(act.e)))
  P$G$attr$i$health = c(t(E$A[1:tf,]))
  # build & save gif
  g = plot.network(P$G,list(fill='health'),list(linetype='active',alpha=1)) +
    labs(title='Day {frame_time}') +
    transition_time(t) +
    ease_aes('linear')
  f = fname(sprintf('gif-n%d-s%d-t%d.gif',N,seed,tf)) # filename
  animate(g,nframes=tf,fps=fps,w=600,h=500,renderer=gifski_renderer(f))
}

# -----------------------------------------------------------------------------
# main

if (sys.nframe() == 0){
  .debug = FALSE
  fit.ptr.engage()
  fit.ptr.kenya()
  fit.ptr.context()
  fit.pdur.kenya()
  fit.durs.plot()
  .debug = TRUE
  E = fit.run.epi()
  fit.plot.epidemic(E)
  fit.plot.tree(E)
  fit.network.gif()
}
