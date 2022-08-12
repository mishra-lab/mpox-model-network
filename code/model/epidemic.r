library('reshape2')

# TODO: rename Z -> x

epi.t = function(t0=1,tf=180){
  t = seq(t0,tf)
}

epi.random.init = function(P,t){
  .Random.seed <<- P$seed.state
  U = list()
  f.sex.e = P$G$attr$e$sex / P$G$attr$g$dur # sex frequency per partnership
  U$e.sex.t = lapply(t,function(tj){ which(runif(P$G$N.e) < f.sex.e) }) # partners had sex per day
  U$u.sex.t = lapply(U$e.sex.t,runif) # random number per sex day
  U$dur.exp.i = P$dur.exp.rfun(P$N) # random durations to onset per-person
  U$dur.inf.i = P$dur.inf.rfun(P$N) # random durations to recovery per-person
  U$dur.iso.i = P$dur.iso.rfun(P$N) # random durations to isolation per-person
  return(U)
}

epi.init.state = function(P){
  S0 = seqn(P$N) # node indices "i"
  I0 = sample(S0,P$N.I0,p=P$G$attr$i$par[S0])
  S0 = setdiff(S0,I0)
  V10 = sample(S0,min(P$N.V0[1],len(S0)))
  S0 = setdiff(S0,V10)
  V20 = sample(S0,min(P$N.V0[2],len(S0)))
  S0 = setdiff(S0,V20)
  X = list()
  X$i = list()
  X$i$S = S0        # i of susceptible
  X$i$E = numeric() # i of exposed
  X$i$I = I0        # i of infected
  X$i$H = numeric() # i of isolating
  X$i$R = numeric() # i of recovered
  X$i$V1 = V10      # i of vaccinated 1 dose
  X$i$V2 = V20      # i of vaccinated 2 dose
  X$dur = list()
  X$dur$exp = numeric()       # duration exposed -> infectious
  X$dur$inf = numeric(P$N.I0) # duration infectious -> {isolating or recovered}
  X$dur$iso = numeric()       # duration infectious -> recovered (among isolating)
  return(X)
}

epi.out.init = function(P,t,X){
  out.t = list()
  out.t$N = list()
  out.t$inc = list()
  out.t$tree = list()
  out.t$Xi = list('t0'=X$i)
  return(out.t)
}

epi.ii.transmit = function(P,U,tj,X){
  # get i of newly infected (exposed)
  ii.sex = P$G$ii.e[U$e.sex.t[[tj]],,drop=FALSE] # partners who had sex today
  b.IZ = ii.sex[,1] %in% X$i$I # IZ partnership (among above)
  b.ZI = ii.sex[,2] %in% X$i$I # ZI partnership
  u.sex = U$u.sex.t[[tj]][c(which(b.IZ),which(b.ZI))] # random number for each partnership
  ii.IZ = rbind(ii.sex[b.IZ,c(1,2)],ii.sex[b.ZI,c(2,1)]) # I-(Z = anything) partnerships
  beta = lookup.map(ii.IZ[,2], X$i, P$beta.health) # beta from health status (among above)
  b.transmit = u.sex < beta
  ii.IE.j = ii.IZ[b.transmit,,drop=FALSE] # infected, exposed
  ii.IE.j = ii.IE.j[!duplicated(ii.IE.j[,2]),,drop=FALSE] # remove double infections
}

epi.run = function(P,t){
  # run the epidemic
  U = epi.random.init(P,t)
  X = epi.init.state(P)
  out.t = epi.out.init(P,t,X)
  for (tj in t){
    # computing transitions
    ii.IE.j = epi.ii.transmit(P,U,tj,X)
    b.EI.j = X$dur$exp > U$dur.exp.i[X$i$E]           # E -> I
    b.IR.j = X$dur$inf > U$dur.inf.i[X$i$I]           # I -> R
    b.IH.j = X$dur$inf > U$dur.iso.i[X$i$I] & !b.IR.j # I -> H
    b.HR.j = X$dur$iso > U$dur.inf.i[X$i$H]           # H -> R
    b.IZ.j = b.IR.j | b.IH.j
    iRj = c(X$i$I[b.IR.j],X$i$H[b.HR.j])
    iHj = X$i$I[b.IH.j]
    iIj = X$i$E[b.EI.j]
    iEj = ii.IE.j[,2]
    iEj.S  = intersect(iEj, X$i$S)
    iEj.V1 = intersect(iEj, X$i$V1)
    iEj.V2 = intersect(iEj, X$i$V2)
    # log state stuff
    out.t$N[[tj]] = lens(X$i)
    out.t$inc[[tj]] = c('S'=len(iEj.S),'V1'=len(iEj.V1),'V2'=len(iEj.V2))
    out.t$tree[[tj]] = ii.IE.j
    # update durations
    X$dur$iso = 1 + c(X$dur$iso[!b.HR.j], X$dur$inf[b.IH.j]) # iso dur continues from inf
    X$dur$inf = 1 + c(X$dur$inf[!b.IZ.j], numeric(len(iIj))) # new inf from zero
    X$dur$exp = 1 + c(X$dur$exp[!b.EI.j], numeric(len(iEj))) # new exp from zero
    # update indices
    X$i$R  = c(X$i$R,iRj)
    X$i$H  = c(X$i$H[!b.HR.j], iHj)
    X$i$I  = c(X$i$I[!b.IZ.j], iIj)
    X$i$E  = c(X$i$E[!b.EI.j], iEj)
    X$i$V2 = setdiff(X$i$V2, iEj.V2)
    X$i$V1 = setdiff(X$i$V1, iEj.V1)
    X$i$S  = setdiff(X$i$S,  iEj.S)
    if (.debug && sum(lens(X$i)) != P$N){ stop('len(X$i) != P$N') }
    if (.debug && lens(X$dur) != lens(X$i[2:4])){ stop('lens(X$dur) != lens(X$i)') }
  }
  out.t$Xi[['tf']] = X$i
  return(epi.results(P,t,out.t))
}

epi.run.s = function(P.s,t,parallel=TRUE){
  # run for multiple seeds, and usually compute the results immediately too
  if (parallel){ lapply.fun = par.lapply } else { lapply.fun = lapply }
  R.s = lapply.fun(P.s,function(P){ R = epi.run(P,t) })
}

epi.results = function(P,t,out.t){
  # collect some results
  R = list()
  P$G = epi.net.attrs(P$G,t,out.t)
  R$P = P
  R$out = epi.output(P,t,out.t)
  R$tree = epi.tree(P,t,out.t)
  return(R)
}

epi.net.attrs = function(G,t,out.t){
  # add some attributes to G after running the model
  G$attr$i$inf.src = factor(G$i %in% out.t$Xi$t0$I,levels=c(F,T),labels=M$inf.src$name)
  G$attr$i$health  = lookup.map(G$i, out.t$Xi$tf)
  return(G)
}

epi.tree = function(P,t,out.t){
  # clean up tree & compute a few properties
  tree = do.call(rbind,lapply(t,function(tj){
    tree.j = out.t$tree[[tj]]
    cbind('par'=tree.j[,1],'chi'=tree.j[,2],'t'=rep(tj,nrow(tree.j)))
  }))
  tree = rbind(cbind('par'=rep(0,P$N.I0),'chi'=out.t$Xi$t0$I,'t'=rep(0,P$N.I0)),tree)
  if (.debug){
    tree.data = recurse.tree(tree[,c('par','chi')],root=0)
    tree = as.data.frame(tree)
    tree = rbind(c(-1,0,NA),tree)
    tree = tree[match(tree.data[1,],tree$chi),]
    tree$gen = tree.data[2,]
    tree$pos = (tree.data[3,]-1) / (max(tree.data[3,])-1)
    tree$n.chi.dir = tree.data[4,]
    tree$n.chi.tot = tree.data[5,]
    tree$dt = tree$t - tree$t[match(tree$par,tree$chi)]
  }
  return(tree)
}

epi.output = function(P,t,out.t){
  # clean-up outputs + compute a few extra
  N       = as.data.frame(do.call(rbind,out.t$N))
  N$all   = P$N
  inc     = as.data.frame(do.call(rbind,out.t$inc))
  inc$all = rowSums(inc)
  prev    = N / P$N
  out  = cbind(
    't' = t,
    setNames(N,   paste0('N.',   names(N))),
    setNames(prev,paste0('prev.',names(prev))),
    setNames(inc, paste0('inc.', names(inc))))
}

epi.output.melt = function(out,P){
  # melt the data in out (usually for plotting)
  N.t = nrow(out)
  out.long = melt(out,id.vars='t')
  out.long = split.col(out.long,'variable',c('var','health'),del=TRUE)
  out.long$seed = P$seed
  return(out.long)
}

epi.output.melt.s = function(R.s){
  # apply epi.output.melt to a list of R.s -- e.g. from epi.run.s
  out.long.s = kw.call(rbind,lapply(R.s,function(R){
    out.long = epi.output.melt(R$out,R$P)
  }))
}