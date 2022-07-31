
.spec = list(
  # some default colors, etc. for plotting
  color = list(
    'city'   = list('A'='#66ccff','B'='#4488dd','C'='#2244bb','bridge'='#cccccc'),
    'health' = list('S'='#ffcc00','I'='#ff0066','R'='#cc00cc','V'='#00cc66'),
    'src'    = list('Import'='#ff0066','Local'='#0099cc')
  ),
  shape = list(
    'src' = list('Import'='square','Local'='circle')
  )
)

.attr.cts = function(G,what,name){
  # get continuous 'what' (node or edge) attribute named 'name'
  if (is.numeric(name)){
    f = name
  } else {
    f = get(paste0(what,'.attr.vec'))(G,name)
  }
}

.attr.map = function(G,what,name,spec){
  # get 'what' (node or edge) attribute named 'name', and remap it to a .spec (e.g. 'color')
  map = .spec[[spec]][[name]]
  if (is.null(map)){
    f = name
  } else {
    a = get(paste0(what,'.attr.vec'))(G,name)
    f = as.character(factor(a,levels=names(map),labels=map))
  }
  return(f)
}

draw.network = function(P,n.color='city',e.color='city',shape='circle',size=NULL){
  # wrapper for plot.igraph, with some default values
  if (is.null(size)){ size = 80/sqrt(len(P$G)) }
  plot(P$G,
    margin       = 0,
    vertex.label = NA,
    edge.color   = .attr.map(P$G,'edge',e.color,'color'),
    vertex.color = .attr.map(P$G,'node',n.color,'color'),
    vertex.shape = .attr.map(P$G,'node',shape,'shape'),
    vertex.size  = .attr.cts(P$G,'node',size),
    layout       = graph_attr(P$G,'layout') # OK if NULL
  )
}

plot.epidemic = function(out.long,y='N',select=list(city='all'),intervals=.9,facet='~',color='health'){
  # plot median for out.long[[y]], maybe after selecting some rows
  # out.long can also be out.long.s (e.g. from rbind), then we add confidence intervals (ci)
  out.long = q3.aggr(y,c('t','health','city'),out.long,intervals=intervals) # compute the ci
  yq = function(q){ paste0(y,'.',q) } # convenience
  ok = rep(TRUE,nrow(out.long)) # selector (initially all)
  for (name in names(select)){
    ok = ok & out.long[[name]] %in% select[[name]] # removing rows ...
  }
  g = ggplot(out.long[ok,],aes_string(x='t')) +
    geom_line(aes_string(y=yq(.5),color=color)) +
    scale_color_manual(values=unlist(.spec$color[[color]]),drop=FALSE) +
    scale_fill_manual( values=unlist(.spec$color[[color]]),drop=FALSE) +
    labs(x='Time (days)') +
    facet_grid(facet) +
    theme_light()
  for (i in intervals){
    ci = q.interval(i)
    g = g + geom_ribbon(aes_string(ymin=yq(ci[1]),ymax=yq(ci[2]),fill=color),alpha=.2)
  }
  return(g)
}
