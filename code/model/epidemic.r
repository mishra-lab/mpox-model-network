library('reshape2')

# TODO: rename U -> R (intuitive); R -> E (); t -> t.vec (transpose conflict)

epi.t = function(t0=1,tf=180){
  t = seq(t0,tf)
}

epi.random.init = function(P,t){
  # U is a list of pre-computed random numbers
  # which allow the exact same events to occur among every individual,
  # unless they are affected by vaccination or other intervention.
  # if we did not pre-compute these numbers, the random number stream
  # would get offset as fewer people get infected (due to interventions)
  # and the exact same events would not occur among other people
  # possibly resulting in "negative intervention impact" due to this extra randomness
  .Random.seed <<- P$seed.state # resume state we left off after def.params
  U = list()
  f.sex.e = P$G$attr$e$sex / P$G$attr$g$dur # sex frequency per partnership (per day)
  U$e.sex.t = lapply(t,function(tj){ which(runif(P$G$N.e) < f.sex.e) }) # partnerships having sex each day
  U$u.sex.t = lapply(U$e.sex.t,runif) # random number per sex-day (used to decide transmission)
  U$dur.EI.i = P$dur.EI.rfun(P$N) # random dur per-person: from exposure to symptom onset (incubation)
  U$dur.IR.i = P$dur.IR.rfun(P$N) # random dur per-person: from onset to recovery (infectious)
  U$dur.IH.i = P$dur.IH.rfun(P$N) # random dur per-person: from onset to isolation (non-isolated)
  return(U)
}

epi.init.state = function(P){
  # state = lists of indices & selected durations
  # e.g. if population looks like c('R','S','I','S'), X$i = list(S=c(2,4),I=c(3),R=c(1))
  # values in X$dur are matched (same order) as corresponding indices in X$i
  # e.g. we couldd have X$dur$I = c(5) for above, if person i=3 became infectious 5 days ago
  S0 = P$G$i # node indices "i"
  I0 = sample(S0,P$N.I0,p=P$G$attr$i$par[S0])
  S0 = setdiff(S0,I0)
  V10 = sample(S0,min(P$N.V0[1],len(S0)))
  S0 = setdiff(S0,V10)
  V20 = sample(S0,min(P$N.V0[2],len(S0)))
  S0 = setdiff(S0,V20)
  X = list()        # state = indices & durations
  X$i = list()      # list of indices in each state
  X$i$S = S0        # i of susceptible
  X$i$E = numeric() # i of exposed
  X$i$I = I0        # i of infected
  X$i$H = numeric() # i of isolating
  X$i$R = numeric() # i of recovered
  X$i$V1 = V10      # i of vaccinated 1 dose
  X$i$V2 = V20      # i of vaccinated 2 dose
  X$dur = list()    # list of durations in selected states
  X$dur$E = numeric()       # duration since exposure before symptom onset (become infectious)
  X$dur$I = numeric(P$N.I0) # duration since onset before {isolating or recovered}
  X$dur$IH = numeric()      # duration since onset before recovered (among isolating)
  return(X)
}

epi.out.init = function(P,t,X){
  # out.t = list of lists; sub-lists indexed by time to allow efficient adding
  # then we can do.call(rbind,out.t${whatever}) to make a data.frame
  out.t = list()
  out.t$N = list()     # number of people in each state -- i.e. lens(X$i)
  out.t$inc = list()   # incidence (count) among S, V1, V2
  out.t$ii.IE = list() # IE (infector -> exposed) pairs
  out.t$Xi = list('t0'=X$i)
  return(out.t)
}

epi.ii.transmit = function(P,U,tj,X){
  # get i of newly infected (exposed)
  ii.sex = P$G$ii.e[U$e.sex.t[[tj]],,drop=FALSE] # partners who had sex today
  b.IA = ii.sex[,1] %in% X$i$I # IA partnership (among above, A = anything)
  b.AI = ii.sex[,2] %in% X$i$I # AI partnership
  u.sex = U$u.sex.t[[tj]][c(which(b.IA),which(b.AI))] # random number for each partnership
  ii.IA = rbind(ii.sex[b.IA,c(1,2)],ii.sex[b.AI,c(2,1)]) # IA partnerships (in that order)
  beta = lookup.map(ii.IA[,2], X$i, P$beta.health) # beta from health status (among above)
  b.transmit = u.sex < beta
  ii.IE.j = ii.IA[b.transmit,,drop=FALSE] # infected, newly exposed
  ii.IE.j = ii.IE.j[!duplicated(ii.IE.j[,2]),,drop=FALSE] # remove double infections
}

epi.run = function(P,t){
  # run the epidemic
  U = epi.random.init(P,t) # pre-compute all (most) random values
  X = epi.init.state(P) # indices of all health states, and durations of some of those
  out.t = epi.out.init(P,t,X) # initialize the output stuff
  for (tj in t){
    # computing transitions
    ii.IE.j = epi.ii.transmit(P,U,tj,X)             # I{S,V1,V2} -> IE
    b.EI.j = X$dur$E  > U$dur.EI.i[X$i$E]           # E -> I (end incubation)
    b.IR.j = X$dur$I  > U$dur.IR.i[X$i$I]           # I -> R (end infectious / recovered)
    b.IH.j = X$dur$I  > U$dur.IH.i[X$i$I] & !b.IR.j # I -> H (begin isolation)
    b.HR.j = X$dur$IH > U$dur.IR.i[X$i$H]           # H -> R (end isolation / recovered)
    b.IA.j = b.IR.j | b.IH.j              # any reason to leave I (H or R)
    iRj = c(X$i$I[b.IR.j],X$i$H[b.HR.j])  # i of newly recovered (from I or H)
    iHj = X$i$I[b.IH.j]                   # i of newly isolating
    iIj = X$i$E[b.EI.j]                   # i of newly infecious
    iEj = ii.IE.j[,2]                     # i of newly exposed
    iEj.S  = intersect(iEj, X$i$S)        # i of newly exposed from S
    iEj.V1 = intersect(iEj, X$i$V1)       # i of newly exposed from V1
    iEj.V2 = intersect(iEj, X$i$V2)       # i of newly exposed from V2
    # log state stuff
    out.t$N[[tj]] = lens(X$i) # how many people in each health state
    out.t$inc[[tj]] = c('S'=len(iEj.S),'V1'=len(iEj.V1),'V2'=len(iEj.V2)) # stratified incidence
    out.t$ii.IE[[tj]] = ii.IE.j # log transmission pairs
    # update durations (remove people who finished this state & append new people)
    X$dur$IH = 1 + c(X$dur$IH[!b.HR.j], X$dur$I[b.IH.j]) # dur$IH continues from dur$I
    X$dur$I  = 1 + c(X$dur$I [!b.IA.j], numeric(len(iIj))) # dur$I restarts from zero
    X$dur$E  = 1 + c(X$dur$E [!b.EI.j], numeric(len(iEj))) # dur$E restarts from zero
    # update indices
    X$i$R  = c(X$i$R,iRj)            # add new recovered
    X$i$H  = c(X$i$H[!b.HR.j], iHj)  # remove recovered, add newly isolating
    X$i$I  = c(X$i$I[!b.IA.j], iIj)  # remove isolating & recovered, add newly infecious
    X$i$E  = c(X$i$E[!b.EI.j], iEj)  # remove newly infectious, add newly exposed
    X$i$V2 = setdiff(X$i$V2, iEj.V2) # remove newly exposed
    X$i$V1 = setdiff(X$i$V1, iEj.V1) # remove newly exposed
    X$i$S  = setdiff(X$i$S,  iEj.S)  # remove newly exposed
    # debug stuff
    if (.debug){
      out.t$Xi[[tj]] = X$i # log full i state every day (expensive)
      if (sum(lens(X$i)) != P$N){ stop('sum(lens(X$i)) != P$N') }
      if (any(lens(X$dur) != lens(X$i[2:4]))){ stop('lens(X$dur) != lens(X$i)') }
    }
  }
  out.t$Xi[['tf']] = X$i
  return(epi.results(P,t,out.t))
}

epi.run.s = function(P.s,t,parallel=TRUE){
  # run for multiple seeds, usually in parallel
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
  if (.debug){
    R$A = dn.array(list('t'=t,'i'=seqn(P$N)),character())
    for (tj in t){ Xij = out.t$Xi[[tj]]; for (h in names(Xij)){ R$A[tj,Xij[[h]]] = h } }
  }
  return(R)
}

epi.net.attrs = function(G,t,out.t){
  # add some attributes to G after running the model
  G$attr$i$inf.src = factor(G$i %in% out.t$Xi$t0$I,levels=c(F,T),labels=M$inf.src$name)
  G$attr$i$health.t0 = lookup.map(G$i, out.t$Xi$t0) # initial health state
  G$attr$i$health.tf = lookup.map(G$i, out.t$Xi$tf) # final health state
  return(G)
}

epi.tree = function(P,t,out.t){
  # clean up tree & compute a few properties
  tree = do.call(rbind,lapply(t,function(tj){
    ii.IE.j = out.t$ii.IE[[tj]]
    cbind('par'=ii.IE.j[,1],'chi'=ii.IE.j[,2],'t'=rep(tj,nrow(ii.IE.j))) # parent (I), child (E), time (t)
  }))
  # add dummy rows for seed cases
  tree = rbind(cbind('par'=rep(0,P$N.I0),'chi'=out.t$Xi$t0$I,'t'=rep(0,P$N.I0)),tree)
  if (.debug){ # recurse.tree can get expensive
    tree.data = recurse.tree(tree[,c('par','chi')],root=0)
    tree = as.data.frame(tree)
    tree = rbind(c(-1,0,NA),tree) # another dummy node (fixes matching)
    tree = tree[match(tree.data[1,],tree$chi),] # re-order to match tree.data
    tree$gen = tree.data[2,] # generation
    tree$pos = (tree.data[3,]-1) / (max(tree.data[3,])-1) # position, rescaled to range [0,1]
    tree$n.chi.dir = tree.data[4,] # number of direct children
    tree$n.chi.tot = tree.data[5,] # number of total direct + indirect children
    tree$dt = tree$t - tree$t[match(tree$par,tree$chi)] # generation time
  }
  return(tree)
}

epi.output = function(P,t,out.t){
  # clean-up outputs + compute a few extra
  N       = as.data.frame(do.call(rbind,out.t$N)) # num people in each health state, each day
  N$all   = P$N
  inc     = as.data.frame(do.call(rbind,out.t$inc)) # num new infections among S,V1,V2 each day
  inc$all = rowSums(inc) # total incidence
  prev    = N / P$N # prevalence
  out  = cbind( # join these data column-wise
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
  out.long$seed = P$seed # in case we rbind data from multiple seeds
  return(out.long)
}

epi.output.melt.s = function(R.s){
  # apply epi.output.melt to a list of R.s -- e.g. from epi.run.s
  out.long.s = kw.call(rbind,lapply(R.s,function(R){
    out.long = epi.output.melt(R$out,R$P)
  }))
}