
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

plot.epidemic = function(out.long,y='N',select=list(city='all'),intervals=.9,facet=NULL,color='health',...){
  # TODO: clean-up intervals using stat_summary / median_hilow (?)
  # plot median for out.long[[y]], maybe after selecting some rows
  # out.long can also be out.long.s (e.g. from rbind), then we add confidence intervals (ci)
  map = list(color=color,...)
  out.long = q3.aggr(y,c('t','health','city'),out.long,intervals=intervals) # compute the ci
  yq = function(q){ paste0(y,'.',q) } # convenience
  g = ggplot(row.select(out.long,select),aes_string(x='t')) +
    geom_line(kw.call(aes_string,map,y=yq(.5))) +
    labs(x='Time (days)') +
    facet_grid(facet) +
    theme_light()
  map.ci = list.update(map,fill=color,color=NULL)
  for (i in intervals){
    ci = q.interval(i)
    g = g + geom_ribbon(kw.call(aes_string,map.ci,ymin=yq(ci[1]),ymax=yq(ci[2])),alpha=.2)
  }
  g = add.meta.scales(g,list.update(map,fill=color))
  return(g)
}
