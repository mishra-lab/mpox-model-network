
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

plot.epidemic = function(out.long,select=list(),conf.int=.9,facet=NULL,color='health',...){
  # plot median for out.long$value, after selecting some rows
  # out.long can also be out.long.s (e.g. from rbind), then we add confidence intervals (ci)
  select = list.update(list(var='N',city='all'),select)
  map = list(color=color,...)
  g = ggplot(row.select(out.long,select),aes_string(x='t',y='value')) +
    stat_summary(geom='line',fun=median,kw.call(aes_string,map)) +
    labs(x='Time (days)') +
    facet_grid(facet) +
    theme_light()
  for (ci in conf.int){
    g = g + stat_summary(geom='ribbon',fun.data=median_hilow,fun.args=list(conf.int=ci),
      kw.call(aes_string,list.update(map,fill=color,color=NULL)),alpha=.2)
  }
  g = add.meta.scales(g,list.update(map,fill=color))
  return(g)
}
