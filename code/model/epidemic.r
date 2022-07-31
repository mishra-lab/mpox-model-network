library('reshape2')

# TODO: negative indices instead of setiff (!)

epi.t = function(t0=1,tf=365){
  seq(t0,tf)
}

epi.init.state = function(P){ # "X"
  i = seq(vcount(P$G)) # node indices
  X = list()
  X$I = sample.strat(i,P$N.I0.city, # i of infected
    strat=node.attr.vec(P$G,'city'),
    weights=node.attr.vec(P$G,'degree'))
  X$Z = sample(seq(P$Z0.max),P$N.I0,rep=TRUE) # days since infection
  X$S = setdiff(i,X$I) # i of susceptible
  X$R = numeric()      # i of recovered
  X$V = numeric()      # i of vaccinated
  return(X)
}

epi.mat.init = function(P,t){
  # large, complete representation: t (rows) x individuals (cols)
  M = array(character(),dim=c(len(t),P$N),dimnames=list('t'=t,'i'=seq(P$N)))
}

epi.mat.update = function(M,X,tj){
  # update M with current state
  M[tj,X$S] = 'S'
  M[tj,X$I] = 'I'
  M[tj,X$R] = 'R'
  M[tj,X$V] = 'V'
  return(M)
}

epi.vaccinate = function(P,X,tj){
  # if within the vaccinating period, select some of X$S to vaccinate
  if ((tj > P$vax.t0) & (tj <= P$vax.t0 + P$vax.dur)) {
    Vj = sample.strat(X$S, P$N.vax.city.day,
      strat=node.attr.vec(P$G,'city')[X$S],
      weights=node.attr.vec(P$G,P$w.vax.attr)[X$S])
  } else {
    Vj = numeric(0)
  }
}

epi.run = function(P,t){
  # run the epidemic
  set.seed(P$seed)
  X = epi.init.state(P)
  M = epi.mat.init(P,t)
  for (tj in t){
    M = epi.mat.update(M,X,tj) # log state
    # if (len(X$I) != len(X$Z)){ stop('len(X$I) != len(X$Z)') } # DEBUG
    # if (len(X$S)+len(X$I)+len(X$R)+len(X$V) != P$N){ stop('len(X) != P$N') } # DEBUG
    Vj  = epi.vaccinate(P,X,tj) # vaccinated TODO: PEP
    nIj = rbinom(len(X$I),degree(P$G)[X$I],P$beta*P$inf.pz[X$Z])
    Ij  = unlist(mapply(sample,adjacent_vertices(P$G,X$I),nIj)) # i of new I
    Ij  = unique(setdiff(Ij,c(X$I,X$R,X$V,Vj))) # TODO: non-perfect vaccine instead
    Zj  = rep(0,len(Ij))        # initialize Z for new I
    bRj = X$Z >= P$inf.dur      # check who has recovered
    X$Z = X$Z + 1               # step days since infection
    X$R = c(X$R,X$I[bRj])       # append new recovered
    X$I = c(X$I[!bRj],Ij)       # remove recovered & append new infected
    X$Z = c(X$Z[!bRj],Zj)       # remove recovered & append new infected
    X$S = setdiff(X$S,c(Ij,Vj)) # remove infected & vaccinated
    X$V = c(X$V,Vj)             # append new vaccinated
  }
  return(M)
}

epi.run.s = function(P.s,t,results=TRUE){
  # run for multiple seeds, and usually compute the results immediately too
  if (results){
    par.lapply(P.s,function(P){ epi.results(P,t,epi.run(P,t)) })
  } else {
    par.lapply(P.s,epi.run(P,t))
  }
}

epi.results = function(P,t,M){
  # collect some results (don't include M, which is large)
  R = list()
  P$G = epi.net.attrs(P$G,M)
  R$P = P
  R$t = t
  R$out = epi.output(P,t,M)
  return(R)
}

epi.net.attrs = function(G,M){
  # add some attributes to G after running the model
  vertex_attr(G,'sqrt.deg') = sqrt(as.numeric(degree(G)))
  vertex_attr(G,'src')      = as.character(factor(M[1,]=='I',labels=c('Local','Import')))
  vertex_attr(G,'health')   = M[length(t),]
  return(G)
}

epi.output = function(P,t,M){
  # sum up the numbers of people in each state over time
  # yields a more efficient representation of M, but loses information
  city = node.attr.vec(P$G,'city')
  out = data.frame('t'=t)
  for (h in c('S','I','R','V')){
    out[[join.str(h,'all')]] = rowSums(M==h)
    for (y in P$lab.city){
      out[[join.str(h,y)]] = rowSums(M[,city==y]==h)
    }
  }
  return(out)
}

epi.output.melt = function(out,P){
  # melt the data in out (usually for plotting)
  N.t = nrow(out)
  out.long = melt(out,id.vars='t')
  out.long = rename.cols(out.long,value='N')
  out.long = split.col(out.long,'variable',c('health','city'),del=TRUE)
  out.long$health = as.character(out.long$health)
  out.long$city   = as.character(out.long$city)
  out.long$seed   = P$seed
  N.city = rep(c(P$N,P$N.city),each=N.t,times=4) # NOTE: 4 from SIRV
  out.long$n.city = out.long$N / N.city # per-city prevalence
  return(out.long)
}

epi.output.melt.s = function(R.s,P.s){
  # apply epi.output.melt to a list of R.s -- e.g. from epi.run.s
  out.long.s = kw.call(rbind,lapply(P.s,function(P){
    out.long = epi.output.melt(R.s[[P$seed]]$out,P)
  }))
}