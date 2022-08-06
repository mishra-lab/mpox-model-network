library('reshape2')

epi.t = function(t0=1,tf=180){
  t = seq(t0,tf)
}

epi.random.init = function(P,t){
  .Random.seed <<- P$seed.state
  dn = list('t'=t,'i'=seqn(P$N))
  N.U = len(t) * P$N
  U = list()
  U$S.E = dn.array(dn,runif(N.U))
  U$E.I = dn.array(dn,runif(N.U))
  U$I.R = dn.array(dn,runif(N.U))
  return(U)
}

epi.init.state = function(P){
  S0 = seqn(P$N) # node indices "i"
  I0 = sample(S0,P$N.I0,p=P$G$attr$i$par[S0])
  S0 = setdiff(S0,I0)
  V10 = sample(S0,P$N.V10,p=P$G$attr$i$par[S0])
  S0 = setdiff(S0,V10)
  V20 = sample(S0,P$N.V20,p=P$G$attr$i$par[S0])
  S0 = setdiff(S0,V20)
  X = list()
  X$S = S0        # i of susceptible
  X$E = numeric() # i of exposed
  X$I = I0        # i of infected
  X$R = numeric() # i of recovered
  X$V1 = V10      # i of vaccinated 1 dose
  X$V2 = V20      # i of vaccinated 2 dose
  return(X)
}

epi.array.init = function(P,t){
  # large, complete representation: t (rows) x individuals (cols)
  A = dn.array(list('t'=t,'i'=seqn(P$N)),character())
}

epi.array.update = function(A,tj,X){
  # update A with current state
  A[tj,X$S]  = 'S'
  A[tj,X$E]  = 'E'
  A[tj,X$I]  = 'I'
  A[tj,X$R]  = 'R'
  A[tj,X$V1] = 'V1'
  A[tj,X$V2] = 'V2'
  return(A)
}

epi.do.expose = function(P,U,tj,X){
  # get i of newly infected
  r = rle(sort(adjacent.i(P$G,X$I)))
  beta = rep(P$beta.health,lengths(X))[match(r$values,unlist(X))]
  Ej = r$values[U$S.E[tj,r$values] < 1-(1-beta/P$net.dur)^r$lengths]
}

epi.do.onset = function(P,U,tj,X){
  # get i of newly infectious / symptomatic (assumed same)
  Ij = X$E[U$E.I[tj,X$E] < 1/P$dur.exp]
}

epi.do.recovery = function(P,U,tj,X){
  # get i of recovered
  Rj = X$I[U$I.R[tj,X$I] < 1/P$dur.inf]
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
    X$R  = c(X$R,Rj)                      # append new recovered
    X$I  = setdiff(c(X$I,Ij),Rj)          # append new infectious & remove recovered
    X$E  = setdiff(c(X$E,Ej),Ij)          # append new exposed & remove infectious
    X$S  = setdiff(X$S,Ej)                # remove exposed
    X$V1 = setdiff(X$V1,Ej)               # remove exposed
    X$V2 = setdiff(X$V2,Ej)               # remove exposed
    if (.debug && sum(lengths(X)) != P$N){ stop('len(X) != P$N') }
  }
  return(epi.results(P,t,A))
}

epi.run.s = function(P.s,t,parallel=TRUE){
  # run for multiple seeds, and usually compute the results immediately too
  if (parallel){ lapply.fun = par.lapply } else { lapply.fun = lapply }
  R.s = lapply.fun(P.s,function(P){ epi.run(P,t) })
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

epi.output.melt.s = function(R.s,P.s){
  # apply epi.output.melt to a list of R.s -- e.g. from epi.run.s
  out.long.s = kw.call(rbind,lapply(P.s,function(P){
    out.long = epi.output.melt(R.s[[P$seed]]$out,P)
  }))
}