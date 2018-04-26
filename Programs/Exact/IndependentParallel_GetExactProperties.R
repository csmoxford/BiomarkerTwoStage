# Functions for independent parallel.

library(EurosarcBayes)

getIndParallelScreenMean = function(p.n,p.p,n1,n,r1,pr.positive,nRep  = 10000){

  prStop.P = sum(dbinom(0:r1,n1,p.p))
  prStop.N = sum(dbinom(0:r1,n1,p.n))


  st = rep(0,nRep)
  for(i in 1:nRep){
    n.p = ifelse(runif(1) < prStop.P,n1,n)
    n.n = ifelse(runif(1) < prStop.N,n1,n)
    marker = rbinom(1000,1,pr.positive)
    st[i] = max(which(cumsum(marker == 1) == n.p)[1],which(cumsum(marker == 0) == n.n)[1])
  }
  st
}

getIndParallelProperties = function(name = "IndParallel",pr.positive, p0,p1, res, n1, n, eta, zeta, prior.a,prior.b){

  res = bayes_binom_one_postprob_nstage(c(n1,n),p0=p0, p1=p1, eta=eta, zeta=zeta, prior.a = prior.a, prior.b = prior.b, round=TRUE)
  r1 = res@failure[1]

  nscreenSample = getIndParallelScreenMean(p0,p0,n1,n,r1,pr.positive)

  dfnames = c("Design","P(+)","P(b-)","P(b+)","interim stop","interim unselected","interim positive","Futility","Success (b- & b+)","Success (b+)","mean SS","99% SS","mean tested")
  df = data.frame(matrix(0,nrow=0,ncol=length(dfnames)), stringsAsFactors = FALSE)
  names(df) = dfnames

  df[1,] = c(name,pr.positive,p0,p0, NA, NA, NA, roundWZero((1-res@alpha)^2,3), roundWZero(res@alpha^2,3), roundWZero((1-res@alpha)*res@alpha,3), roundWZero(2*res@exp.p0,1) , 2*n, roundWZero(mean(nscreenSample),1))

  nscreenSample = getIndParallelScreenMean(p0,p1,n1,n,r1,pr.positive)

  df[2,]  = c(name,pr.positive,p0,p1, NA, NA, NA, roundWZero((1-res@alpha)*(1-res@power),3), roundWZero(res@alpha*res@power,3), roundWZero((1-res@alpha)*res@power,3), roundWZero(res@exp.p0 + res@exp.p1,1), 2*n, roundWZero(mean(nscreenSample),1))

  nscreenSample = getIndParallelScreenMean(p1,p1,n1,n,r1,pr.positive)

  df[3,]  = c(name,pr.positive,p1,p1, NA, NA, NA, roundWZero((1-res@power)^2,3), roundWZero(res@power^2,3), roundWZero((1-res@power)*res@power,3),roundWZero(2*res@exp.p1,1), 2*n, roundWZero(mean(nscreenSample),1))

  df
}
