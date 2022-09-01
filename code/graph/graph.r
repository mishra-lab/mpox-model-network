source('utils/all.r')

# TODO: add comments per-function
# TODO: change notation i -> n ("node" more cononical)

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
  if (is.null(deg.i)){ deg.i = degrees.from.edges(i,ii.e) }
  G = list()
  G$N.i = len(i)
  G$N.e = nrow(ii.e)
  G$i = i
  G$e = seqn(G$N.e)
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
  if (N.i == 0){ return(matrix(nrow=0,ncol=2)) }
  if (shuffle){ i = sample(i) }
  ii.e = edges.low.high(cbind(i[1:(N.i/2)],i[(N.i/2+1):N.i]))
}

edges.group.odds = function(i,g,or.gg,shuffle=TRUE){
  if (shuffle){ ord = order(sample(i)); i = i[ord]; g = g[ord] }
  i.g = split(i,g) # split i by g
  N.g = lens(i.g)  # count i by g
  u.g = 1:len(N.g) # unique levels in g
  N.gg.0 = outer(N.g,N.g)/sum(N.g)/2 # proportional mixing
  N.gg = iter.prop.fit(or.gg*N.gg.0,N.g/2) # applying odds (preferences)
  for (g in u.g){ N.gg[g,] = round.sum(N.gg[g,],N.g[g]/2) } # sum-preserving rounding
  i.gg = lapply(u.g,function(g){ # assign i to each gg case
    # e.g. i.gg[[1]][[2]] is the i in group 1 having edges with group 2
    split(i.g[[g]],factor(rep(u.g,N.gg[g,]),u.g))
  })
  gg.m = combn(u.g,2)
  ii.e = do.call(rbind,c( # build edges
    lapply(u.g,function(g){ # for g1 = g2
      ii.e.g = edges.unloop(edges.random(i.gg[[g]][[g]],shuffle=FALSE))
    }),
    lapply(1:ncol(gg.m),function(m){ # for g1 != g2
      g1 = gg.m[1,m]; g2 = gg.m[2,m];
      ii.e.gg = edges.unloop(edges.random(c(i.gg[[g1]][[g2]],i.gg[[g2]][[g1]]),shuffle=FALSE))
    })
  ))
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
# tree stuff

.tree.tips <<- 0 # TODO: pass around internally?
tree.recurse = function(ii,root=0,gen=0,pos=NULL){
  # assuming ii represents a tree, walk the tree (in given ordder) & return matrix with columns:
  # index (ordered by tree search), generation, position, n direct children, n total children
  if (is.null(pos)){ pos = 'child.range' }
  b.root = ii[,1]==root
  i.childs = ii[b.root,2]
  if (any(b.root)){
    mat.childs = matrix(nrow=5,unlist(lapply(i.childs,function(i.child){
      tree.recurse(ii=ii[!b.root,,drop=FALSE],root=i.child,gen=gen+1,pos=pos)
    })))
    root.pos = switch(pos,
      'tip.mean' = mean(mat.childs[3,]),
      'tip.range' = mean(range(mat.childs[3,])),
      'child.mean' = mean(mat.childs[3,mat.childs[1,] %in% i.childs]),
      'child.range' = mean(range(mat.childs[3,mat.childs[1,] %in% i.childs])))
    mat.root.childs = matrix(nrow=5,c(root,gen,root.pos,sum(b.root),ncol(mat.childs),mat.childs))
  } else {
    .tree.tips <<- .tree.tips + 1
    mat.root = matrix(nrow=5,c(root,gen,.tree.tips,0,0))
  }
}

tree.pc.map = function(tree,cols,par='par',chi='chi'){
  # most row data in a tree corresponds to the child (in our implementation, anyways)
  # here we lookup the parent of each child, and copy some data from those rows
  # onto the same row as the child, creating new columns with suffix .par
  # usually for plotting
  pc.map = match(tree[[par]],tree[[chi]])
  for (col in cols){
    tree[[paste0(col,'.',par)]] = tree[[col]][pc.map]
  }
  return(tree)
}

# ==================================================================================================
# layout + plotting

graph.layout.random = function(G){
  # layout nodes randomly (uniform) within unit circle, no consideration of edges
  r = sqrt(runif(G$N.i))
  t = runif(G$N.i) * 2 * pi
  pos = matrix(c(r*cos(t),r*sin(t)),ncol=2)
}

graph.layout.fr = function(G,niter=500,temp.pwr=.5,temp.tau=10,dc.pwr=1,gif.pos=NULL){
  # layout nodes using Fruchterman-Reingold algorithm (attract-repulse)
  temp = G$N.i^temp.pwr # initial speed
  rtemp = (1-log(2)*temp.tau/niter) # speed decay
  dc = G$N.i^dc.pwr # controls global spread
  iie = as.factor(c(G$ii.e)) # pre-compute for speed
  pos = graph.layout.random(G) # random initial position
  b.tri = lower.tri(matrix(0,G$N.i,G$N.i),diag=TRUE) # pre-compute
  tri.0 = function(x){ x[b.tri] = 0; x }
  for (k in 1:niter){
    # if (k %% 10 == 0){ G$attr$g$layout = pos; print(plot.graph(G)) } # DEBUG
    if (is.list(gif.pos)){ gif.pos[[k]] = pos } # DEBUG: gif
    # distances - TODO: this can blow up, is there a faster/less-mem way?
    d1 = outer(pos[,1],pos[,1],'-') # dx
    d2 = outer(pos[,2],pos[,2],'-') # dy
    dd = d1^2 + d2^2
    # repulsive forces
    ddr = (dc - dd * sqrt(dd)) / (dd * dc)
    d1r = tri.0(d1) * ddr
    d2r = tri.0(d2) * ddr
    d1r.sum = rowSums(d1r,na.rm=TRUE) - colSums(d1r,na.rm=TRUE)
    d2r.sum = rowSums(d2r,na.rm=TRUE) - colSums(d2r,na.rm=TRUE)
    # attractive forces
    d1e = d1[G$ii.e]
    d2e = d2[G$ii.e]
    dde = dd[G$ii.e]
    d1e.sum = sapply(split(c(-d1e,d1e)*dde,iie),sum)
    d2e.sum = sapply(split(c(-d2e,d2e)*dde,iie),sum)
    # update pos & temp
    dpos = cbind(d1r.sum + d1e.sum, d2r.sum + d2e.sum)
    pos = pos + dpos * temp / sqrt(dpos[,1]^2 + dpos[,2]^2)
    temp = temp * rtemp
  }
  if (is.list(gif.pos)){ return(gif.pos) } # return gif.pos if using
  pos = pos - rep(colMeans(pos),each=G$N.i) # return centred
}

plot.graph = function(G,i.aes=list(),e.aes=list()){
  if (is.null(G$attr$g$layout)){ G$attr$g$layout = graph.layout.fr(G) }
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
  g = ggplot(data=i.data,aes(x=X,y=Y)) +
    kw.call(geom_curve,data=e.data,e.aes[!b.e],
      map=kw.call(aes_string,e.aes[b.e],x='X1',y='Y1',xend='X2',yend='Y2')) +
    kw.call(geom_point,data=i.data,i.aes[!b.i],
      kw.call(aes_string,i.aes[b.i])) +
    guides(fill=guide_legend(override.aes=list(shape=21))) +
    coord_equal() +
    theme_void()
}

# ==================================================================================================
# casting to igraph if needed

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
