suppressPackageStartupMessages({library('igraph')})

fitness.game = function(n,deg.mean,fitness,adj.power){
  G = static.fitness.game(
    no.of.edges = round(n * deg.mean),
    fitness.out = fitness
  )
  p.add = degree(G)^adj.power / (1+degree(G)) # empiric
  G = connect.iso(G,p.add=p.add)
}

gen.net = function(net.type,N,net.args,node.attrs,edge.attrs){
  # generate a network 
  net.fun = switch(net.type,
    'er'  = random.graph.game,
    'ba'  = barabasi.game,
    'ff'  = forest.fire.game,
    'fit' = fitness.game)
  G = kw.call(net.fun,net.args,n=N)
  vertex_attr(G,'degree') = degree(G)
  for (name in names(node.attrs)){
    vertex_attr(G,name) = node.attrs[name]
  }
  for (name in names(edge.attrs)){
    edge_attr(G,name) = edge.attrs[name]
  }
  return(G)
}

gen.multi.city = function(P){ # TODO: port to graph/graph.r
  # generate a network from multiple cities
  city.attr = lapply(P$lab.city,function(city){ list('city'=city) })
  G.city.list = mapply(gen.net,P$net.type,P$N.city,P$net.args,city.attr,city.attr)
  G = do.call(graph.disjoint.union,G.city.list)
  graph_attr(G,'layout') = layout_nicely(G) # DEBUG: pre-compute consistent network layout
  e.scale.bridge = round(2 * ecount(G) * P$e.bridge / sum(P$mix.bridge)) # scale mix.bridge
  G = bridge.nets(G,'city',P$mix.bridge*e.scale.bridge)
  # print(mean(E(G)$city=='bridge')) # DEBUG: should ~= P$e.brideg
  return(G)
}

bridge.nets = function(G,g.attr,mix.bridge){ # TODO: port to graph/graph.r
  # if G contains multiple sub-nets, defined by edge attribute g.attr (e.g. 'city')
  # this function breaks some edges within each group & reconnects them per mix.bridge
  # mix.bridge must be symmetric, zero diagonal, and rowSums divisible by 2 exactly
  e.swap = rowSums(mix.bridge)
  if (!isSymmetric(mix.bridge) | any(diag(mix.bridge)!=0) | any((e.swap %% 2)>0)){ stop('Invalid mix.bridge') }
  E.mat = get.edgelist(G) # edge matrix, size: ecount(G) x 2
  E.g   = edge.attr.vec(G,g.attr) # vec of attributes by which to group & bridge
  g.vec = unique(E.g) # groups
  n.g   = len(g.vec)
  g.seq = seqn(n.g)
  E.mat.del = lapply(g.seq,function(g){ # randomly select edges to delete within each group
    E.mat.g = E.mat[E.g==g.vec[g],]
    matrix(E.mat.g[sample(nrow(E.mat.g),size=e.swap[g]/2),],ncol=2)
  })
  E.mat.add = matrix(nrow=0,ncol=2)
  e.del = lapply(E.mat.del,function(E.mat.g){ c(E.mat.g) }) # vector of available nodes per group
  for (g1 in g.seq){
    for (g2 in g.seq){
      if ((g2 < g1) & (mix.bridge[g1,g2] > 0)){ # nonzero upper triangle elements
        e.slice = seqn(mix.bridge[g1,g2]) # indices of nodes within e.del
        E.mat.add = rbind(E.mat.add,cbind(e.del[[g1]][e.slice],e.del[[g2]][e.slice])) # add new edge
        e.del[[g1]] = e.del[[g1]][-e.slice] # remove available nodes per group
        e.del[[g2]] = e.del[[g2]][-e.slice] # remove available nodes per group
  }}}
  G = del.edges(G,kw.call(rbind,E.mat.del))
  G = kw.call(add.edges,setNames(list('bridge'),g.attr),G,E.mat.add)
  return(G)
}

connect.iso = function(G,p.add=NULL){ # TODO: port to graph/graph.r (?)
  # re-connect isolated nodes to non-isolated nodes (with weights p.add)
  N = vcount(G)
  if (N == 0){ return(G) }
  i = seqn(N)
  b.0 = degree(G) == 0
  i.0 = i[b.0]
  i.add = sample(i[!b.0],len(i.0),p=p.add[!b.0])
  G = add.edges(G,c(rbind(i.0,i.add)))
}

node.attr.vec = function(G,name){
  # return as vector instead of list
  a = unname(unlist(vertex_attr(G,name)))
}

edge.attr.vec = function(G,name){
  # return as vector instead of list
  a = unname(unlist(edge_attr(G,name)))
}

del.edges = function(G,E.mat){
  # see https://github.com/igraph/rigraph/issues/550
  G = delete_edges(G,apply(E.mat,1,paste0,collapse='|'))
}

add.edges = function(G,E.mat,...){
  # see https://github.com/igraph/rigraph/issues/550
  G = add_edges(G,c(t(E.mat)),...)
}
