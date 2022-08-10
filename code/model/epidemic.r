library('reshape2')

# TODO: rename X$N -> X$i

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
  X$i$R = numeric() # i of recovered
  X$i$V1 = V10      # i of vaccinated 1 dose
  X$i$V2 = V20      # i of vaccinated 2 dose
  X$dur = numeric(P$N) # duration in current state
  return(X)
}

epi.array.init = function(P,t){
  # large, complete representation: t (rows) x individuals (cols)
  A = dn.array(list('t'=t,'i'=seqn(P$N)),character())
}

epi.array.update = function(A,tj,X){
  # update A with current state
  A[tj,X$i$S]  = 'S'
  A[tj,X$i$E]  = 'E'
  A[tj,X$i$I]  = 'I'
  A[tj,X$i$R]  = 'R'
  A[tj,X$i$V1] = 'V1'
  A[tj,X$i$V2] = 'V2'
  return(A)
}

epi.do.expose = function(P,U,tj,X){
  # get i of newly infected (exposed)
  ii.sex = matrix(P$G$ii.e[U$e.sex.t[[tj]],],ncol=2) # partners who had sex today
  b.IZ = ii.sex[,1] %in% X$i$I # IZ partnership (among above)
  b.ZI = ii.sex[,2] %in% X$i$I # ZI partnership
  u.sex = U$u.sex.t[[tj]][c(which(b.IZ),which(b.ZI))] # random number for each partnership
  i.Z = c(ii.sex[b.IZ,2],ii.sex[b.ZI,1]) # i of partners of I (may also be I)
  beta = NA * i.Z # initialize beta (among above) - susceptibility
  for (h in names(X$i)){ beta[i.Z %in% X$i[[h]]] = P$beta.health[h] } # beta from health status
  Ej = unique(i.Z[u.sex < beta]) # infected (exposed)
}

epi.do.onset = function(P,U,tj,X){
  # get i of newly infectious / symptomatic (assumed same)
  Ij = X$i$E[X$dur[X$i$E] > U$dur.exp.i[X$i$E]]
}

epi.do.recovery = function(P,U,tj,X){
  # get i of recovered (or isolating)
  Rj = X$i$I[X$dur[X$i$I] > U$dur.inf.i[X$i$I] | X$dur[X$i$I] > U$dur.iso.i[X$i$I]]
}

epi.run = function(P,t){
  # run the epidemic
  U = epi.random.init(P,t)
  X = epi.init.state(P)
  A = epi.array.init(P,t)
  for (tj in t){
    A = epi.array.update(A,tj,X) # log state
    # computing transitions
    Ej = epi.do.expose(P,U,tj,X)
    Ij = epi.do.onset(P,U,tj,X)
    Rj = epi.do.recovery(P,U,tj,X)
    # applying transitions
    X$i$R  = c(X$i$R,Rj)                      # append new recovered
    X$i$I  = setdiff(c(X$i$I,Ij),Rj)          # append new infectious & remove recovered
    X$i$E  = setdiff(c(X$i$E,Ej),Ij)          # append new exposed & remove infectious
    X$i$S  = setdiff(X$i$S,Ej)                # remove exposed
    X$i$V1 = setdiff(X$i$V1,Ej)               # remove exposed
    X$i$V2 = setdiff(X$i$V2,Ej)               # remove exposed
    # update durations
    X$dur[c(Ej,Ij,Rj)] = 0
    X$dur = X$dur + 1
    if (.debug && sum(lengths(X$i)) != P$N){ stop('len(X$i) != P$N') }
  }
  return(epi.results(P,t,A))
}

epi.run.s = function(P.s,t,parallel=TRUE){
  # run for multiple seeds, and usually compute the results immediately too
  if (parallel){ lapply.fun = par.lapply } else { lapply.fun = lapply }
  R.s = lapply.fun(P.s,function(P){ R = epi.run(P,t) })
}

epi.results = function(P,t,A){
  # collect some results (usually don't include A, which is large)
  R = list()
  P$G = epi.net.attrs(P$G,A,t)
  R$P = P
  R$t = t
  R$out = epi.output(P,t,A)
  if (.debug){ R$A = A }
  return(R)
}

epi.net.attrs = function(G,A,t){
  # add some attributes to G after running the model
  G$attr$i$inf.src = factor(A[1,]=='I',levels=c(F,T),labels=M$inf.src$name)
  G$attr$i$health  = A[length(t),]
  return(G)
}

epi.output = function(P,t,A){
  # computs some pre-determined outputs from A, usually stratified by health
  # yields a more efficient representation of A, but loses information
  # TODO: maybe save A (compressed somehow?) to avoid information loss + compute outputs on the fly
  out = data.frame('t'=t)
  # function compute outputs overall (as list)
  sum.fun = function(num,den=1,...){
    out[[join.str(...)]] = rowSums(num) / { if (is.matrix(den)) rowSums(den) else den }
    return(out)
  }
  # absolute (N), prevalence (prev)
  A1 = matrix(1,nrow=nrow(A),ncol=ncol(A))
  out = sum.fun(A!='',1,'N','all')
  for (h in M$health$name){
    out = sum.fun(A==h,1,'N',h)
    out = sum.fun(A==h,A1,'prev',h)
  }
  # absolute incidence (inc)
  h.sus = c('S','V1','V2')
  A.1 = A[(2:len(t))-1,]
  A.E = A[(2:len(t)),] == 'E'
  out = sum.fun(rbind(A.1 %in% h.sus & A.E,NA),1,'inc','all')
  for (h in h.sus){
    out = sum.fun(rbind(A.1==h & A.E,NA),1,'inc',h)
  }
  return(out)
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