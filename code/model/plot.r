
aggr.quantile = function(out.long,x,y,map,facet=NULL,ci=.9){
  qci = c(.5,(1-ci)/2,1-(1-ci)/2)
  if (!is.null(facet)){ map = c(map,strsplit(facet,'\\s?[~+]\\s?')[[1]]) }
  strat = map[map %in% colnames(out.long)]
  f = paste(y,'~',paste(c(x,strat),collapse=' + '))
  out.long.q = do.call(data.frame,aggregate(formula(f),out.long,quantile,qci,names=F))
}

aes.quantile = function(x,y,map){
  do.call(aes_string,c(list(x=x,y=paste0(y,'.1'),ymin=paste0(y,'.2'),ymax=paste0(y,'.3')),map))
}

plot.clean = function(g){
  g = g + theme_light() +
    theme(strip.background=element_rect(fill='gray90'),strip.text=element_text(color='black'))
}

add.meta.scales = function(g,map){
  for (name in names(map)){
    value = map[[name]]
    if (value %in% names(M)){
      Mi = M[[value]]
      scale = get(paste0('scale_',name,'_manual'))
      g = g + scale(values=Mi[[name]],labels=Mi$label,name=Mi$title)
    }
  }
  return(g)
}

plot.network = function(G,i.aes=list(),e.aes=list()){
  g = plot.graph(G,i.aes=i.aes,e.aes=e.aes)
  g = add.meta.scales(g,i.aes)
  g = add.meta.scales(g,e.aes)
}

plot.tree = function(tree,...,t.root=NA,pc.map=TRUE){
  # plot a transmission tree vs time
  if (pc.map){ # if faceting, need to pre-compute pc.map for each facet
    tree$t[1] = t.root # use negative value to connect I0 in the past
    tree = tree.pc.map(tree,c('t','pos'))
  }
  args = list(...)
  map = list.update(list(fill='black',color='white',shape=21),args)
  b.aes = map %in% names(tree)
  g = ggplot(tree) +
    geom_segment(aes_string(y='t.par',x='pos.par',xend='pos',yend='t'),alpha=.5) +
    geom_point(aes_string(color=map[b.aes]$fill,x='pos',y='t'),size=0) +
    kw.call(geom_point,map=kw.call(aes_string,map[b.aes],x='pos',y='t'),map[!b.aes]) +
    scale_x_continuous(labels=NULL,breaks=NULL) +
    labs(y='Time (days)',x='')
  g = plot.clean(g)
}

add.tree.margin = function(g,fill=TRUE){
  # add marginal to tree, but must do this last as ggMarginal creates a gtable
  g = g + theme(legend.position='top')
  g = ggExtra::ggMarginal(g,type='hist',margins='y',
    lwd=.25,color='white',alpha=.8,groupFill=!is.null(fill),
    yparams=list(binwidth=1))
}

plot.durs = function(rfuns,n=1e5,xmax=30,dx=1){
  x    = seq(0,xmax,dx)
  data = do.call(rbind,lapply(names(rfuns),function(name){
    y = rfuns[[name]](n)
    pdf = hist(y,c(x-dx/2,Inf),plot=FALSE)$density * dx
    var = factor(rep(1:2,each=len(x)),1:2,c('Duration','Survival Function'))
    value = c(pdf,1-cumsum(pdf))
    data.frame(x=x,value=value,var=var,name=name)
  }))
  g = ggplot(data,aes(x=x,y=value,color=name,fill=name)) +
    geom_bar(stat='identity',alpha=.2,width=dx,color=NA,show.legend=FALSE) +
    geom_step(position=position_nudge(x=-dx/2),show.legend=FALSE) +
    facet_grid('var ~ name',scales='free_y') +
    labs(x='Time (Days)',y='Distribution Value')
  g = plot.clean(g)
}

plot.G.distr = function(P.s,gie,attr,select=list(),...,cuts=NULL,cum=FALSE){
  # plot distribution(s) of attributes in P$G (or P.s)
  # note: no geom by default for flexibility; we just collect the histogram data
  if ('G' %in% names(P.s)){ P.s = list(P.s) }
  f = formula(paste('n','~',paste(c(attr,'seed',...),collapse=' + ')))
  data.raw = do.call(rbind,lapply(P.s,function(P){
    cbind('seed'=P$seed,'n'=1,as.data.frame(P$G$attr[[gie]]))
  }))
  data = aggregate(f,row.select(data.raw,select),sum)
  if (!is.null(cuts)){ data[[attr]] = int.cut(data[[attr]],cuts) }
  if (cum){
    ids = lapply(list('seed',...),function(v){ data[[v]] })
    data$n = do.call(ave,c(list(x=data$n,FUN=cum1),ids)) }
  g = ggplot(data,aes_string(x=attr,y='n',...))
  g = add.meta.scales(g,list(...))
  g = plot.clean(g)
}

plot.epidemic = function(out.long,select=list(),conf.int=.9,facet=NULL,color='health',scales=NULL,...){
  # plot median for out.long$value, after selecting some rows
  # out.long can also be out.long.s (e.g. from rbind), then we add ribbon confidence interval
  out.long = row.select(out.long,list.update(list(var='N'),select))
  map = list(color=color,fill=color,...)
  out.long.q = aggr.quantile(out.long,x='t',y='value',map=map,facet=facet,ci=conf.int)
  g = ggplot(out.long.q,aes.quantile(x='t',y='value',map=map)) +
    geom_ribbon(color=NA,alpha=.2) +
    geom_line() +
    labs(x='Time (days)') +
    facet_grid(facet,scales=scales)
  g = add.meta.scales(g,list.update(map,fill=color))
  g = plot.clean(g)
  return(g)
}
