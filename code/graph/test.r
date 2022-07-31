source('graph/graph.r')

def.arg.test = function(N,seed=NULL){
  arg = list()
  arg$seed = seed
  arg$N = N
  arg$deg.shape = 0.255
  arg$deg.rate = 0.032
  arg$deg.shift = 0.913
  arg$main.frac = .2
  arg$main.p.deg.power = -1
  arg$main.m.shape = 5
  arg$main.m.rate = .2
  arg$casu.m.shape = .3
  arg$casu.m.rate = .3
  return(arg)
}

make.graph.test = function(arg,debug=FALSE){
  set.seed(arg$seed)
  i = seqn(arg$N)
  # sample degrees
  deg.i = round(rgamma(arg$N,shape=arg$deg.shape,arg$deg.rate) + arg$deg.shift)
  deg.i = degrees.balanced(deg.i)
  # generate main partners
  i.main = sample(i,round(arg$main.frac*len(i)),p=deg.i^arg$main.p.deg.power)
  ii.e.main = edges.random(i.main,shuffle=FALSE)
  r.e.main = round(rgamma(nrow(ii.e.main),shape=arg$main.m.shape,rate=arg$main.m.rate))
  ii.e.main.r = edges.repeated(ii.e.main,r.e.main)
  # generate casual partnerss
  deg.i.casu = deg.i
  deg.i.casu[i.main] = deg.i.casu[i.main] - 1
  ii.e.casu = edges.unloop(edges.from.degrees(i,deg.i.casu))
  r.e.casu = round(1+rgamma(nrow(ii.e.casu),shape=arg$casu.m.shape,rate=arg$casu.m.rate))
  ii.e.casu.r = edges.repeated(ii.e.casu,r.e.casu)
  # all contacts
  ii.e = rbind(ii.e.main.r,ii.e.casu.r)
  # attributes
  g.attr = list()
  i.attr = list()
  e.attr = list()
  if (debug){
    i.attr$par.p6m = deg.i
    i.attr$sex.p6m = degrees.from.edges(i,ii.e)
    i.attr$main = factor(i %in% i.main,c(T,F),c('Yes','No'))
    e.attr$type = factor(c(rep('main',nrow(ii.e.main.r)),rep('casu',nrow(ii.e.casu.r))))
    kr.e = index.repeated.edges(ii.e)
    e.attr$k.e = kr.e[,1]
    e.attr$r.e = kr.e[,2]
  }
  # graph object
  G = graph.obj(ii.e=ii.e,i=i,deg.i=deg.i,g.attr=g.attr,i.attr=i.attr,e.attr=e.attr)
  if (debug){ G$attr$g$layout = graph.layout(G) }
  return(G)
}

test.speed = function(){
  file.out = file.path('.tmp',paste0('graph-benchmark-',Sys.Date(),'.out'))
  out = microbenchmark::microbenchmark(
    {    make.graph.test(def.arg.test(100),debug=TRUE) },
    {   make.graph.test(def.arg.test(100),debug=FALSE) },
    {  make.graph.test(def.arg.test(1000),debug=FALSE) },
    { make.graph.test(def.arg.test(10000),debug=FALSE) },
  times=100)
  print(out); sink(file.out); print(out); sink()
}

test.plot = function(N=500,S=1){
  adj = list(
    fill  = scale_fill_viridis(option='inferno',trans='log10',limits=c(1,180)),
    color = scale_fill_viridis(option='viridis',trans='log10',limits=c(1,180)),
    size = scale_size_continuous(range=c(.5,3),breaks=c(1,3,10,30)))
  for (s in seqn(S)){
    G = make.graph.test(def.arg.test(N,seed=s),debug=TRUE)
    print(plot.graph(G,list(fill='par.p6m')) + adj$fill + labs(fill='P6M\nPartners'))
    print(plot.graph(G,list(fill='sex.p6m')) + adj$fill + adj$size + labs(fill='P6M\nSex'))
  }
}

test.speed()
test.plot()