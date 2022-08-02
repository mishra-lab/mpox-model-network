library('reshape2')

# TODO: negative indices instead of setiff, faster ?

epi.t = function(t0=1,tf=180){
  t = seq(t0,tf)
}

epi.init.state = function(P){ # "X"
  S0 = seqn(P$N) # node indices "i"
  I0 = sample.strat(S0,P$N.I0.city,
    strat   = P$G$attr$i$city,
    weights = P$G$attr$i$par.p6m)
  S0 = setdiff(S0,I0)
  V10 = epi.do.vaccinate(P,list(S=S0,V1=S0),0)[[1]]
  S0 = setdiff(S0,V10)
  V20 = epi.do.vaccinate(P,list(S=S0,V1=S0),0)[[2]]
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
  A = array(character(0),dim=c(len(t),P$N),dimnames=list('t'=t,'i'=seqn(P$N)))
}

epi.array.update = function(A,X,tj){
  # update A with current state
  A[tj,X$S]  = 'S'
  A[tj,X$E]  = 'E'
  A[tj,X$I]  = 'I'
  A[tj,X$R]  = 'R'
  A[tj,X$V1] = 'V1'
  A[tj,X$V2] = 'V2'
  return(A)
}

epi.do.vaccinate = function(P,X,tj){
  # get i of newly vaccinated
  Vj = list(numeric(),numeric()) # dose 1, dose 2
  for (P.vax in P$vax.params.phase){ # for each vaccination phase
    b.tj = P.vax$t == tj # when is tj within P.vax$t
    if (any(b.tj)){
      i = switch(P.vax$dose,'1'=X$S,'2'=X$V1) # dose -> sampling from S or V1
      w = switch(is.null(P.vax$w.attr),T=NULL,F=P$G$attr$i[[P.vax$w.attr]][i]) # weights
      N.tj.city = P.vax$N.day.city[b.tj,] # number of doses today per city
      Vj.phase = sample.strat(i,N.tj.city,strat=P$G$attr$i$city[i],weights=w) # sample
      Vj[[P.vax$dose]] = c(Vj[[P.vax$dose]],Vj.phase) # append vax from this phase 
    }
  }
  Vj = lapply(Vj,unique) # remove any duplicates
}

epi.do.breakthrough = function(P,X){
  # get i of vaccinated who could experience breakthrough - TODO: can this be done (faster) via p?
  VSj = c(sample.i(X$V1,1-P$vax.eff.dose[1]),sample.i(X$V2,1-P$vax.eff.dose[2]))
}

epi.do.expose = function(P,X,Sj){
  # get i of newly infected
  Ej = sample.i(unlist(adjacent.i(P$G,X$I)),P$beta/P$net.dur)
  Ej = unique(intersect(Ej,Sj)) # unique exposed who are susceptible
}

epi.do.infectious = function(P,X){
  # get i of infectious
  Ij = sample.i(X$E,1/P$dur.exp)
}

epi.do.recovery = function(P,X){
  # get i of recovered
  Rj = sample.i(X$I,1/P$dur.inf)
}

epi.run = function(P,t){
  # run the epidemic
  .Random.seed <<- P$seed.state
  X = epi.init.state(P)
  A = epi.array.init(P,t)
  for (tj in t){
    A = epi.array.update(A,X,tj) # log state
    # computing transitions
    Vj = epi.do.vaccinate(P,X,tj) # TODO: PEP
    V1j = Vj[[1]]; V2j = Vj[[2]];
    Sj = c(X$S,epi.do.breakthrough(P,X))
    Ej = epi.do.expose(P,X,Sj)
    Ij = epi.do.infectious(P,X)
    Rj = epi.do.recovery(P,X)
    # applying transitions
    X$R  = c(X$R,Rj)                      # append new recovered
    X$I  = setdiff(c(X$I,Ij),Rj)          # append new infectious & remove recovered
    X$E  = setdiff(c(X$E,Ej),Ij)          # append new exposed & remove infectious
    X$V1 = setdiff(c(X$V1,V1j),c(Ej,V2j)) # append new dose-1 & remove exposed, dose-2
    X$V2 = setdiff(c(X$V2,V2j),Ej)        # append new dose-2 & remove exposed
    X$S  = setdiff(X$S,c(Ej,V1j))         # remove exposed, dose-1
    if (.debug && sum(sapply(X,len)) != P$N){ stop('len(X) != P$N') }
  }
  return(A)
}

epi.run.s = function(P.s,t,results=TRUE,parallel=TRUE){
  # run for multiple seeds, and usually compute the results immediately too
  if (parallel){ lapply.fun = par.lapply } else { lapply.fun = lapply }
  if (results){
    lapply.fun(P.s,function(P){ epi.results(P,t,epi.run(P,t)) })
  } else {
    lapply.fun(P.s,function(P){ epi.run(P,t) })
  }
}

epi.results = function(P,t,A){
  # collect some results (don't include A, which is large)
  R = list()
  P$G = epi.net.attrs(P$G,A)
  R$P = P
  R$t = t
  R$out = epi.output(P,t,A)
  return(R)
}

epi.net.attrs = function(G,A){
  # add some attributes to G after running the model
  G$attr$i$inf.src = factor(A[1,]=='I',levels=c(F,T),labels=M$inf.src$name)
  G$attr$i$health  = A[length(t),]
  return(G)
}

epi.output = function(P,t,A){
  # computs some pre-determined outputs from A, usually stratified by health & city
  # yields a more efficient representation of A, but loses information
  # TODO: maybe save A (compressed somehow?) to avoid information loss + compute outputs on the fly
  out = data.frame('t'=t)
  # function compute outputs overall + by city (as list)
  sum.city.fun = function(num,den=1,...){
    out[[join.str(...,'all')]] = rowSums(num) / { if (is.matrix(den)) rowSums(den) else den }
    for (city in M$city$name){
      b.city = P$G$attr$i$city==city
      out[[join.str(...,city)]] = rowSums(num[,b.city]) / { if (is.matrix(den)) rowSums(den[,b.city]) else den }
    }
    return(out)
  }
  # absolute (N), prevalence (prev)
  A1 = matrix(1,nrow=nrow(A),ncol=ncol(A))
  out = sum.city.fun(A!='',1,'N','all')
  for (h in M$health$name){
    out = sum.city.fun(A==h,1,'N',h)
    out = sum.city.fun(A==h,A1,'prev',h)
  }
  # absolute incidence (inc)
  h.sus = c('S','V1','V2')
  A.1 = A[(2:len(t))-1,]
  A.E = A[(2:len(t)),] == 'E'
  out = sum.city.fun(rbind(A.1 %in% h.sus & A.E,NA),1,'inc','all')
  for (h in h.sus){
    out = sum.city.fun(rbind(A.1==h & A.E,NA),1,'inc',h)
  }
  return(out)
}

epi.output.melt = function(out,P){
  # melt the data in out (usually for plotting)
  N.t = nrow(out)
  out.long = melt(out,id.vars='t')
  out.long = split.col(out.long,'variable',c('var','health','city'),del=TRUE)
  out.long$seed = P$seed
  return(out.long)
}

epi.output.melt.s = function(R.s,P.s){
  # apply epi.output.melt to a list of R.s -- e.g. from epi.run.s
  out.long.s = kw.call(rbind,lapply(P.s,function(P){
    out.long = epi.output.melt(R.s[[P$seed]]$out,P)
  }))
}