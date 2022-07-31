source('utils/all.r')

# key notation:
# G: a graph object
# i: indices of nodes (people) in the graph
# e: indices of edges (contacts) in the graph
# N.{y}: total number of index {y}
# deg.i: degree of node i
# ii.e: edge list (matrix with ncol=2) with rows like [a,b] for edge a--b
# r.e: number of times edge e is repeated
# k.e: index of edge e among the repeats (1,2,...,r.e)
# G$attr.{y}: attribute vectors assumed to match the index {y}
# {x}.{y}: thing {x} stratified or indexed by {y}, usually

# ==================================================================================================
# pseudo-class

graph.obj = function(ii.e,i=NULL,deg.i=NULL,g.attr=NULL,i.attr=NULL,e.attr=NULL){
  if (is.null(i)){ i = sort(unique(c(ii.e))) }
  if (is.null(deg.i)){ deg.i = degree.from.edges(i,ii.e) }
  G = list()
  G$N.i = len(i)
  G$N.e = nrow(ii.e)
  G$i = i
  G$ii.e = ii.e
  G$deg.i = deg.i
  G$attr = list()
  G$attr$g = g.attr
  G$attr$i = i.attr
  G$attr$e = e.attr
  return(G)
}

# ==================================================================================================
# graph generation tools

degrees.balanced = function(deg.i,by=2){
  rem = sum(deg.i) %% by
  if (rem == 0){ return(deg.i) }
  i.fix = sample.int(len(deg.i),rem)
  deg.i[i.fix] = deg.i[i.fix] + 1
  return(deg.i)
}

degrees.from.edges = function(i,ii.e){
  deg.i = c(unname(table(factor(c(ii.e),i))))
}

edges.random = function(i,shuffle=TRUE){
  N.i = len(i)
  if (shuffle){ i = sample(i) }
  ii.e = edges.low.high(cbind(i[1:(N.i/2)],i[(N.i/2+1):N.i]))
}

edges.from.degrees = function(i,deg.i){
  ii.e = edges.random(rep(i,times=deg.i))
}

edges.unloop = function(ii.e){
  e.loops = (ii.e[,1] == ii.e[,2])
  n.loops = sum(e.loops)
  if (n.loops == 0){ return(ii.e) }
  if (n.loops == 1){ return(ii.e[!e.loops,]) }
  if (len(unique(ii.e[e.loops,1]))==1){ return(ii.e[!e.loops,]) }
  ii.e[e.loops,2] = sample(ii.e[e.loops,2])
  return(edges.unloop(edges.low.high(ii.e)))
}

edges.repeated = function(ii.e,r.e){
  if (nrow(ii.e)==0){ return(ii.e) }
  ii.e = ii.e[rep(1:nrow(ii.e),times=r.e),]
}

index.repeated.edges = function(ii.e){
  ord.e = edges.sort.order(ii.e)
  ii.e.ord = ii.e[ord.e,]
  e.hash = ii.e.ord[,1] * 1e6 + ii.e.ord[,2]
  R = rle(e.hash)
  k.e.ord = unname(do.call(c,lapply(R$lengths,seqn)))
  r.e.ord = rep(R$lengths,times=R$lengths)
  kr.e = cbind(k.e.ord[order(ord.e)],r.e.ord[order(ord.e)])
}

edges.low.high = function(ii.e){
  return(cbind(pmin(ii.e[,1],ii.e[,2]),pmax(ii.e[,1],ii.e[,2])))
}

edges.sort.order = function(ii.e){
  return(order(ii.e[,1],ii.e[,2]))
}

# ==================================================================================================
# graph usage tools

adjacent.i = function(G,i){
  # TODO: should this just accept ii.e directly?
  i.adj = c(G$ii.e[G$ii.e[,1] %in% i,2], G$ii.e[G$ii.e[,2] %in% i,1])
}

# ==================================================================================================
# igraph stuff + plotting

.igraph = function(G){
  i.deg.0 = which(G$deg.i==0)
  if (len(i.deg.0) == 0){
    iG = igraph::graph_from_edgelist(G$ii.e,dir=FALSE)
  } else {
    ii.e.tmp = cbind(i.deg.0,1)
    ii.e = rbind(ii.e.tmp,G$ii.e)
    iG = igraph::graph_from_edgelist(ii.e,dir=FALSE)
    iG = igraph::delete_edges(iG,1:len(i.deg.0))
  }
}

graph.layout = function(G,...){
  # TODO: keep disconnected groups a little closer?
  layout = igraph::layout_with_fr(.igraph(G),...)
}

plot.graph = function(G,i.aes=list(),e.aes=list()){
  if (is.null(G$attr$g$layout)){ G$attr$g$layout = graph.layout(G) }
  # defaults
  i.aes = list.update(list(size=25/G$N.i^.4,shape=21),xu=i.aes)
  e.aes = list.update(list(curvature=0,alpha=.1),xu=e.aes)
  # required data
  i.data = setNames(data.frame(G$attr$g$layout),c('X','Y'))
  e.data = data.frame(
    X1 = i.data$X[G$ii.e[,1]],
    X2 = i.data$X[G$ii.e[,2]],
    Y1 = i.data$Y[G$ii.e[,1]],
    Y2 = i.data$Y[G$ii.e[,2]])
  # add G$attr to data, split aes into attr vs fixed values
  if (len(G$attr$e)){ e.data = cbind(e.data,G$attr$e) }
  if (len(G$attr$i)){ i.data = cbind(i.data,G$attr$i) }
  b.i = i.aes %in% names(G$attr$i)
  b.e = e.aes %in% names(G$attr$e)
  g = ggplot() +
    kw.call(geom_curve,data=e.data,e.aes[!b.e],
      map=kw.call(aes_string,e.aes[b.e],x='X1',y='Y1',xend='X2',yend='Y2')) +
    kw.call(geom_point,data=i.data,i.aes[!b.i],
      kw.call(aes_string,i.aes[b.i],x='X',y='Y')) +
    guides(fill=guide_legend(override.aes=list(shape=21))) +
    coord_equal() +
    theme_void()
}
