source('utils/all.r')
source('graph/graph.r')

matrix.string = function(x,fmt='%d'){
  x.str = paste(collapse='\n',apply(x,1,function(xi){
    paste(collapse=' ',lapply(xi,sprintf,fmt=fmt))
  }))
}

get.args = function(N.g,or.gg,k.g){
  #   N.g: number of nodes per group
  # or.gg: odds matrix of edges forming between groups
  #   k.g: number of edges per node in each group
  m = len(N.g)
  if (missing(k.g)){ k.g = rep(1,m) }
  k.i = rep(k.g,N.g)
  g.i = rep(1:m,N.g)
  g.e = rep(1:m,N.g*k.g)
  i.i = seqn(sum(N.g))
  i.e = rep(i.i,k.i)
  data = list(
    i.i = i.i,
    i.e = i.e,
    g.i = factor(g.i),
    g.e = factor(g.e),
    or.gg = matrix(or.gg,m,m))
}

test.edges.group.odds = function(){
  args = list(
    get.args(c(80,20),0),
    get.args(c(50,50),0),
    get.args(c(50,50),diag(2)),
    get.args(c(50,50),5*diag(2)),
    get.args(c(50,50),5*(1-diag(2))),
    get.args(c(40,30,20),0),
    get.args(c(40,30,20),c(5,0,0, 0,0,3, 0,3,0)),
    get.args(c(50,50),0,c(1,2)),
    get.args(c(50,50),3*diag(2),c(1,2)),
    get.args(c(50,50),0,c(1,3)),
    get.args(c(50,50),3*diag(2),c(1,3))
  )
  for (a in args){
    ii.e = do.call(edges.group.odds,list(i=a$i.e,g=a$g.e,or.gg=a$or.gg))
    G = graph.obj(ii.e,i=a$i.i,i.attr=list(g=a$g.i))
    go = plot.graph(G,list(fill='g')) +
      scale_fill_viridis(option='inferno',discrete=TRUE,begin=.2,end=.9) +
      ggtitle(matrix.string(a$or.gg))
    print(go)
  }
}

test.edges.group.odds()