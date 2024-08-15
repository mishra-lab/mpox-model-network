source('utils/all.r')
source('model/meta.r')
source('model/plot.r')

data.fig.save = function(var,...){
  fig.save(root.path('out','fig','.tmp','fit','data',paste0('kenya-',var)),...)
}

load.data = function(var){
  X = read.csv(root.path('data','.private','kenya',sprintf('hivst_%s.csv',var)))
}

main.ptr = function(plot=FALSE){
  X = load.data('ptr_x')
  X$cp = ave(X$p,X$county,X$dt,FUN=cumsum) # cum prob
  if (plot){
    g = plot.ptr(X); data.fig.save('ptr',w=7,h=5) }
  return(X)
}

plot.ptr = function(X){
  var.ftr = c('p'='Proportion','cp'='Cumulative proportion')
  X$dts = sprintf('Recall period: %d days',X$dt)
  i = X$county != 'all'
  X$x[i & X$x==6] = NA
  X$x[X$dt==7 & X$x>8] = NA
  X = melt(X,measure=names(var.ftr))
  X$var = factor(X$variable,names(var.ftr),var.ftr)
  g = ggplot(X,aes(y=value,x=x,fill=county,color=county)) +
    facet_grid('var~dts',scales='free') +
    geom_line() + geom_point(size=1) +
    # geom_bar(stat='identity',position=position_dodge(preserve='single'),color=NA) +
    labs(x='Number of sexual partners',y='Value') +
    ylim(c(0,NA))
  g = plot.clean(add.meta.scales(g,list(color='county',fill='county')))
}

main.pdur = function(){
  X = melt(load.data('pdur_x'),id=c('county','p1m','k'),var='m',value.name='n')
  X$dur = as.numeric(gsub('m','',X$m)); X$m = NULL
  X$p1m.c = int.cut(X$p1m,c(0,2,5))
  X$src = 'Data'
  X$p  = ave(X$n,X$county,X$p1m,X$k,FUN=sum1)
  X$cp = ave(X$n,X$county,X$p1m,X$k,FUN=cum1)
  X = X[order(X$p1m,X$k,X$dur),]
  durs = unique(X$dur)*30
  efun = function(par){
    cp = do.call(pweibull,c(list(q=durs),as.list(par)))
    e = sum(X$n*abs(X$cp-cp)^2) }
  opt = optim(c(1,360),efun,method='L-BFGS-B',lower=1e-6,upper=1e6)
  print(opt)
  X.mod = data.frame(county=NA,p1m=NA,p1m.c=NA,k=NA,n=NA,src='Model',dur=durs/30,
    p  = do.call(dweibull,c(list(x=durs),as.list(opt$par))),
    cp = do.call(pweibull,c(list(q=durs),as.list(opt$par))))
  g = plot.pdur(rbind(X,X.mod)); data.fig.save('pdur',w=6,h=4)
}

args.pdur = list(shape=.597,scale=233)

plot.pdur = function(X){
  g = ggplot(X,aes(x=dur,y=cp,color=p1m.c,fill=p1m.c,linetype=src)) +
    stat_summary(geom='ribbon',fun.min=min,fun.max=max,alpha=.3,color=NA) +
    stat_summary(geom='line',fun=median) +
    scale_x_continuous(trans='log10') +
    scale_color_viridis_d(na.value='black') +
    scale_fill_viridis_d(na.value='black') +
    labs(x='Duration (months)',y='Cumulative Proportion',linetype='',
      color='P1M\nPartners',fill='P1M\nPartners') +
    ylim(c(0,NA))
  g = plot.clean(g)
}

# main.ptr(plot=TRUE)
# main.pdur()
