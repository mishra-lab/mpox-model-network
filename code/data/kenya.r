source('utils/all.r')
source('model/meta.r')
source('model/plot.r')

data.fig.save = function(var,...){
  fig.save(root.path('out','fig','.tmp','fit',paste0('kenya-',var)),...)
}

load.data = function(var){
  X = read.csv(root.path('data','.private','kenya',sprintf('hivst_%s.csv',var)))
}

main.ptr = function(plot=FALSE){
  X = load.data('ptr_x')
  X$cp = ave(X$p,X$county,X$dt,FUN=cumsum) # cum prob
  if (plot){
    g = plot.ptr(X); data.fig.save('ptr',w=6,h=5) }
  return(X)
}

plot.ptr = function(X,xmin=6,xmax=20,tform='identity'){
  X$dts = sprintf('Recall period: %d days',X$dt)
  i = X$county != 'all'
  Y = cbind(X[i & X$x >= xmin-1,],xmin=xmin-1:0,xmax=c(xmin-1,xmax))
  X$x[i & X$x==xmin] = NA
  g = ggplot(X,aes(y=cp,x=x,fill=county,color=county)) +
    facet_grid('dts') +
    geom_line() +
    geom_ribbon(data=Y,aes(xmin=xmin,xmax=xmax),color=NA,alpha=.3) +
    scale_x_continuous(trans=tform) +
    labs(x='Number of sexual partners',y='Cumulative proportion') +
    ylim(c(-.001,1.001))
  g = plot.clean(add.meta.scales(g,list(color='county',fill='county')))
}

# main.ptr(plot=TRUE)
