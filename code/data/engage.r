source('utils/all.r')
source('model/meta.r')
source('model/plot.r')

data.fig.save = function(var,...){
  fig.save(root.path('out','fig','.tmp','fit',paste0('engage-',var)),...)
}

load.data = function(var){
  X = read.csv(root.path('data','.private','engage',sprintf('engage_%s.csv',var)))
  X = rename.cols(X,reltn_type='type',rel_status='status',p_rel.rds='p')
}

clean.q = function(X,types){
  X = remove.cols(X,'n_indv','mean','sd')
  X = melt(X,id=c('city','type','n_rel'),var='p') # melt
  X$p = as.numeric(gsub('q_','',X$p)) # clean p labels
  return(X[X$type %in% types,]) # only some ptr types
}

main.q = function(var,label,vscale=1,vmax=1){
  X = clean.q(load.data(paste0(var,'_q')),c('casual','main'))
  # fit parametric
  P = do.call(rbind,lapply(unique(X$type),function(type){
    i = X$p < .99 & X$type == type # rows to consider (all 3 cities)
    ps = X$p[i]     # quantile probs  (data)
    vx = X$value[i] # quantile values (data)
    efun = function(par){ # error function
      pp = pweibull(vx,shape=par[1],scale=vscale*par[2]) # probs for vx
      e = sum(abs(ps-pp)^2) } # error = (data probs - param probs)^2
    opt = optim(c(.5,1),efun,method='L-BFGS-B',lower=1e-6,upper=1e6)
    print(opt$par*c(1,vscale))
    data.frame(p=ps,type=type,value=qweibull(ps,opt$par[1],vscale*opt$par[2]))
  }))
  # plot
  g = ggplot(X,aes(x=value/vscale,y=p)) +
    facet_grid('type') +
    geom_line(aes(color=city)) +
    geom_line(data=P,linetype='22') +
    labs(x=label,y='Cumulative proportion') +
    coord_cartesian(xlim=c(0,vmax))
  g = plot.clean(add.meta.scales(g,list(color='city')))
  data.fig.save(var,w=6,h=5)
}
p.uncl.excl = 0.0 # assign uncl: 50% to excl & 50% to open
stats = list(
  excl = 'main-exclusive',
  open = 'main-open',
  uncl = 'main-unclear',
  noma = 'no-main')

which.stat = function(X){ i = lapply(stats,function(stat){ which(X$status==stat) }) }

main.stat = function(){
  X = load.data('stat_p')
  X = aggregate(p~status,X,mean)
  i = which.stat(X)
  # split-up uncl:
  X[i$excl,'p'] = X[i$excl,'p'] + X[i$uncl,'p'] * (p.uncl.excl)
  X[i$open,'p'] = X[i$open,'p'] + X[i$uncl,'p'] * (1-p.uncl.excl)
  X = X[-i$uncl,]
  lapply(i,function(i){X[i,'p']}) # return props as list
}

main.p6m = function(p,plot=FALSE){
  X = load.data('p6m_x')
  X = rename.cols(remove.cols(X,'n','n.rds'),p.rds='p')
  i = which.stat(X)
  f = function(xi,xu,pi,pu){ (xi*pi+xu*pu)/(pi+pu) } # weighted avg with * and uncl
  X[i$excl,'p'] = mapply(f,X[i$excl,'p'],X[i$uncl,'p'],pi=p$excl,pu=p$uncl*(p.uncl.excl))
  X[i$open,'p'] = mapply(f,X[i$open,'p'],X[i$uncl,'p'],pi=p$open,pu=p$uncl*(1-p.uncl.excl))
  X = X[-i$uncl,]
  X = rbind(X,cbind(city='avg',aggregate(p~x+status,X,mean))) # proportion
  X$cp = ave(X$p,X$status,X$city,FUN=cumsum) # cum prop
  if (plot){
    g = plot.p6m(X,tform='identity'); data.fig.save('p6m-raw',w=6,h=5)
    g = plot.p6m(X,tform='log10');    data.fig.save('p6m-log',w=6,h=5) }
  X$status = factor(X$status,stats,names(stats))
  return(X)
}

plot.p6m = function(X,xmax=100,tform='identity'){
  g = ggplot(X[X$x<xmax,],aes(y=cp,x=x,color=city)) +
    facet_grid('status') +
    geom_line(aes(color=city)) +
    scale_x_continuous(trans=tform) +
    labs(x='Number of sexual partners (past 6 months)',y='Cumulative proportion')
  g = plot.clean(add.meta.scales(g,list(color='city')))
}

# main.q('pdur',label='Duration (years)',vscale=365,vmax=30)
# main.q('fsex',label='Frequency (per day)',vscale=1,vmax=1)
# p = main.stat()
# x = main.p6m(p,plot=TRUE)
