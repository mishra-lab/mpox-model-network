source('utils/all.r')
source('graph/graph.r')

matrix.string = function(x,fmt='%d'){
  x.str = paste(collapse='\n',apply(x,1,function(xi){
    paste(collapse=' ',lapply(xi,sprintf,fmt=fmt))
  }))
}

get.args = function(N.g,or){
  m = len(N.g)
  data = list(
    i = seqn(sum(N.g)),
    g = factor(c(rep(1:len(N.g),N.g))),
    or.gg = matrix(or,m,m)
  )
}

test.edges.group.odds = function(){
  args.s = list(
    get.args(c(50,50), 0),
    get.args(c(50,50), diag(2)),
    get.args(c(50,50), 5*diag(2)),
    get.args(c(50,50), 5*(1-diag(2))),
    get.args(c(80,20), 0),
    get.args(c(40,30,20),0),
    get.args(c(40,30,20),c(5,0,0, 0,0,3, 0,3,0))
  )
  for (args in args.s){
    ii.e = do.call(edges.group.odds,args)
    G = graph.obj(ii.e,i.attr=list(g=args$g))
    go = plot.graph(G,list(fill='g')) +
      scale_fill_viridis(option='inferno',discrete=TRUE,begin=.2,end=.9) +
      ggtitle(matrix.string(args$or.gg))
    print(go)
  }
}

test.edges.group.odds()