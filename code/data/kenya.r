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

# main.ptr(plot=TRUE)
