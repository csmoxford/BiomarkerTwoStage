


meanScreened = function(
  pr.pos = 0.5,
  stopInt,contPos,contBoth,
  n1Neg,n1Pos,
  nNeg,nPos) {
  
  
  ssInt = 0
  prPos = dnbinom(0:(n1Neg-1),n1Pos,pr.pos)
  for(i in 1:length(prPos)) {
    ssInt = ssInt + prPos[i]*(n1Pos + i-1 + (n1Neg - i + 1)/(1 - pr.pos))
  }
  prNeg = dnbinom(0:(n1Pos-1),n1Neg,1 - pr.pos)
  for(i in 1:length(prNeg)) {
    ssInt = ssInt + prNeg[i]*(n1Neg + i-1 + (n1Pos - i + 1)/(pr.pos))
  }
  
  
  ssPosStageTwo = (nPos-n1Pos)/ pr.pos
  
  
  ssBothStageTwo = 0
  prPos = dnbinom(0:(nNeg-n1Neg-1),nPos-n1Pos,pr.pos )
  for(i in 1:length(prPos)) {
    ssBothStageTwo = ssBothStageTwo + prPos[i] * (nPos - n1Pos + i - 1 + (nNeg - n1Neg - i + 1)/(1 - pr.pos))
  }
  prNeg = dnbinom(0:(nPos-n1Pos-1),nNeg-n1Neg,1 - pr.pos )
  for(i in 1:length(prNeg)) {
    ssBothStageTwo = ssBothStageTwo + prNeg[i] * (nNeg - n1Neg + i - 1 + (nPos - n1Pos - i + 1)/(pr.pos))
  }
  
  
  return(ssInt + ssPosStageTwo*contPos + ssBothStageTwo*contBoth)
  
}


ParallelTwoStageBestFirst = function(
  pNeg,pPos, pr.pos = 0.5,
  n1Neg = 17,n1Pos = 17,
  nNeg = 34,nPos = 34,
  r1NegFail = 1,r1PosFail = 1,
  rNegFail = 6,rPosFail = 6) {
  
  stopInt = 0
  contIntPosFail = 0
  contIntPosSuccess = 0
  contIntBothFail = 0
  contIntBothSPos = 0
  contIntBothSBoth = 0
  
  probPos1 = dbinom(0:nPos,n1Pos,pPos)
  probNeg1 = dbinom(0:nNeg,n1Neg,pNeg)
  probPos2 = dbinom(0:(nPos-n1Pos),nPos-n1Pos,pPos)
  probNeg2 = dbinom(0:(nNeg-n1Neg),nNeg-n1Neg,pNeg)
  for(r1Neg in 0:n1Neg) {
    for(r1Pos in 0:n1Pos) {
      for(rPos in 0:(nPos-n1Pos)) {
        for(rNeg in 0:(nNeg-n1Neg)) {
          prob = probPos1[r1Pos+1] * probNeg1[r1Neg+1] * probPos2[rPos+1] * probNeg2[rNeg+1]
          
          # parallel best first
          if(r1Pos <= r1PosFail) { # stop interim all
            stopInt = stopInt + prob
          } else if(r1Neg <= r1NegFail) { # continue +ve only
            if(r1Pos+rPos <= rPosFail) {
              contIntPosFail = contIntPosFail + prob
            } else {
              contIntPosSuccess = contIntPosSuccess + prob 
            }
          } else { # continue all
            if(r1Pos+rPos <= rPosFail) {
              contIntBothFail = contIntBothFail + prob
            } else if(r1Neg + rNeg <= rNegFail) {
              contIntBothSPos = contIntBothSPos + prob
            } else {
              contIntBothSBoth = contIntBothSBoth + prob
            }
          }
          
        }
      }
    }
  }
  
  # sample size
  ss = (n1Neg + n1Pos) + (nPos-n1Pos)*(1-stopInt) + (nNeg-n1Neg)*(contIntBothFail + contIntBothSPos + contIntBothSBoth)
  # h0
  h0 = stopInt + contIntPosFail + contIntBothFail
  # h1
  h1 = contIntPosSuccess + contIntBothSPos
  # h2
  h2 = contIntBothSBoth

  muScreen = meanScreened(pr.pos, stopInt, (contIntPosFail + contIntPosSuccess), (contIntBothFail + contIntBothSPos + contIntBothSBoth), n1Neg, n1Pos, nNeg, nPos)
  
  
  return(list(h0=h0, h1=h1, h2=h2,ss=ss, muScreen = muScreen, intStop = stopInt, intPos = (contIntPosSuccess + contIntPosFail),intBoth = (contIntBothFail + contIntBothSPos + contIntBothSBoth)))
}


ParallelTwoStageWorstFirst = function(
  pNeg,pPos, pr.pos = 0.5,
  n1Neg = 17,n1Pos = 17,
  nNeg = 34,nPos = 34,
  r1NegFail = 1,r1PosFail = 1,
  rNegFail = 6,rPosFail = 6) {
  
  stopInt = 0
  contIntPosFail = 0
  contIntPosSuccess = 0
  contIntBothFail = 0
  contIntBothSPos = 0
  contIntBothSBoth = 0
  
  probPos1 = dbinom(0:nPos,n1Pos,pPos)
  probNeg1 = dbinom(0:nNeg,n1Neg,pNeg)
  probPos2 = dbinom(0:(nPos-n1Pos),nPos-n1Pos,pPos)
  probNeg2 = dbinom(0:(nNeg-n1Neg),nNeg-n1Neg,pNeg)
  for(r1Neg in 0:n1Neg) {
    for(r1Pos in 0:n1Pos) {
      for(rPos in 0:(nPos-n1Pos)) {
        for(rNeg in 0:(nNeg-n1Neg)) {
          prob = probPos1[r1Pos+1] * probNeg1[r1Neg+1] * probPos2[rPos+1] * probNeg2[rNeg+1]
          
          # parallel worst first
          if(r1Neg > r1NegFail) { # continue all
            if(r1Neg + rNeg > rNegFail) {
              contIntBothSBoth = contIntBothSBoth + prob
            } else if(r1Pos+rPos > rPosFail) {
              contIntBothSPos = contIntBothSPos + prob
            } else {
              contIntBothFail = contIntBothFail + prob
            }
          } else if(r1Pos > r1PosFail) { # continue +ve only
            if(r1Pos+rPos > rPosFail) {
              contIntPosSuccess = contIntPosSuccess + prob 
            } else {
              contIntPosFail = contIntPosFail + prob
            }
          } else { # stop
            stopInt = stopInt + prob
          }
        }
      }
    }
  }
  
  # sample size
  ss = (n1Neg + n1Pos) + (nPos-n1Pos)*(1-stopInt) + (nNeg-n1Neg)*(contIntBothFail + contIntBothSPos + contIntBothSBoth)
  # h0
  h0 = stopInt + contIntPosFail + contIntBothFail
  # h1
  h1 = contIntPosSuccess + contIntBothSPos
  # h2
  h2 = contIntBothSBoth
  
  
  muScreen = meanScreened(pr.pos, stopInt, (contIntPosFail + contIntPosSuccess), (contIntBothFail + contIntBothSPos + contIntBothSBoth), n1Neg, n1Pos, nNeg, nPos)
  
  
  return(list(h0=h0, h1=h1, h2=h2,ss=ss, muScreen = muScreen, intStop = stopInt, intPos = (contIntPosSuccess + contIntPosFail),intBoth = (contIntBothFail + contIntBothSPos + contIntBothSBoth)))
}



ParallelTwoStageBalanced = function(
  pNeg,pPos, pr.pos = 0.5,
  n1Neg = 17,n1Pos = 17,
  nNeg = 34,nPos = 34,
  rIntAllFail = 5, rFinalAllFail = 11,
  r1NegFail = 1,r1PosFail = 1,
  rNegFail = 6,rPosFail = 6) {
  
  stopInt = 0
  contIntPosFail = 0
  contIntPosSuccess = 0
  contIntBothFail = 0
  contIntBothSPos = 0
  contIntBothSBoth = 0
  
  probPos1 = dbinom(0:nPos,n1Pos,pPos)
  probNeg1 = dbinom(0:nNeg,n1Neg,pNeg)
  probPos2 = dbinom(0:(nPos-n1Pos),nPos-n1Pos,pPos)
  probNeg2 = dbinom(0:(nNeg-n1Neg),nNeg-n1Neg,pNeg)
  for(r1Neg in 0:n1Neg) {
    for(r1Pos in 0:n1Pos) {
      for(rPos in 0:(nPos-n1Pos)) {
        for(rNeg in 0:(nNeg-n1Neg)) {
          prob = probPos1[r1Pos+1] * probNeg1[r1Neg+1] * probPos2[rPos+1] * probNeg2[rNeg+1]
          if(r1Neg + r1Pos > rIntAllFail) { # interim unselected test
            if(r1Neg > r1NegFail) { # continue all
              #################################################
              # Final analysis with all
              if(r1Neg + r1Pos + rNeg + rPos > rFinalAllFail) {
                if(r1Neg + rNeg > rNegFail) {
                  contIntBothSBoth = contIntBothSBoth + prob
                } else if (r1Pos + rPos > rPosFail) {
                  contIntBothSPos = contIntBothSPos + prob
                } else {
                  contIntBothFail = contIntBothFail + prob
                }
              } else if(r1Pos + rPos > rPosFail) {
                contIntBothSPos = contIntBothSPos + prob
              } else {
                contIntBothFail = contIntBothFail + prob
              }
              
              
              #################################################
            } else if(r1Pos > r1PosFail) { # continue +ve only
              if(r1Pos+rPos > rPosFail) {
                contIntPosSuccess = contIntPosSuccess + prob 
              } else {
                contIntPosFail = contIntPosFail + prob
              }
            } else { # stop
              stopInt = stopInt + prob
            }
            
          } else if(r1Pos > r1PosFail) { # continue +ve only
            if(r1Pos+rPos > rPosFail) {
              contIntPosSuccess = contIntPosSuccess + prob 
            } else {
              contIntPosFail = contIntPosFail + prob
            }
          } else { # stop
            stopInt = stopInt + prob
          }
        }
      }
    }
  }
  
  # sample size
  ss = (n1Neg + n1Pos) + (nPos-n1Pos)*(1-stopInt) + (nNeg-n1Neg)*(contIntBothFail + contIntBothSPos + contIntBothSBoth)
  # h0
  h0 = stopInt + contIntPosFail + contIntBothFail
  # h1
  h1 = contIntPosSuccess + contIntBothSPos
  # h2
  h2 = contIntBothSBoth
  
  muScreen = meanScreened(pr.pos, stopInt, (contIntPosFail + contIntPosSuccess), (contIntBothFail + contIntBothSPos + contIntBothSBoth), n1Neg, n1Pos, nNeg, nPos)
  
  
  return(list(h0=h0, h1=h1, h2=h2,ss=ss, muScreen = muScreen, intStop = stopInt, intPos = (contIntPosSuccess + contIntPosFail),intBoth = (contIntBothFail + contIntBothSPos + contIntBothSBoth)))
}


# approach is one of balanced, worstFirst, bestFirst
ParallelTwoStageProperties = function(
  approach="balanced",p0 = 0.1,p1 = 0.3, pr.pos = 0.5,
  param=list(n1Neg = 17,n1Pos = 17,
  nNeg = 34,nPos = 34,
  rIntAllFail = 5, rFinalAllFail = 11,
  r1NegFail = 1,r1PosFail = 1,
  rNegFail = 6,rPosFail = 6)) {
  
  if(!approach %in% c("balanced","worstFirst","bestFirst")) {
    stop("approach must be one of balanced, worstFirst, bestFirst")
  }
  
  if(approach == "balanced") {
    res = with(param,{
      nn = ParallelTwoStageBalanced(pNeg = p0,pPos = p0, pr.pos = pr.pos,
                                    n1Neg = n1Neg,n1Pos = n1Pos,
                                    nNeg = nNeg,nPos = nPos,
                                    rIntAllFail = rIntAllFail, rFinalAllFail = rFinalAllFail,
                                    r1NegFail = r1NegFail,r1PosFail = r1PosFail,
                                    rNegFail = rNegFail,rPosFail = rNegFail)
      
      np = ParallelTwoStageBalanced(pNeg = p0,pPos = p1, pr.pos = pr.pos,
                                    n1Neg = n1Neg,n1Pos = n1Pos,
                                    nNeg = nNeg,nPos = nPos,
                                    rIntAllFail = rIntAllFail, rFinalAllFail = rFinalAllFail,
                                    r1NegFail = r1NegFail,r1PosFail = r1PosFail,
                                    rNegFail = rNegFail,rPosFail = rNegFail)
      
      pp = ParallelTwoStageBalanced(pNeg = p1,pPos = p1, pr.pos = pr.pos,
                                    n1Neg = n1Neg,n1Pos = n1Pos,
                                    nNeg = nNeg,nPos = nPos,
                                    rIntAllFail = rIntAllFail, rFinalAllFail = rFinalAllFail,
                                    r1NegFail = r1NegFail,r1PosFail = r1PosFail,
                                    rNegFail = rNegFail,rPosFail = rNegFail)
      
      return(list(nn=nn,np=np,pp=pp))
    })
    
  } else if(approach == "worstFirst") {
    res = with(param,{
      nn = ParallelTwoStageWorstFirst(pNeg = p0,pPos = p0, pr.pos = pr.pos,
                                    n1Neg = n1Neg,n1Pos = n1Pos,
                                    nNeg = nNeg,nPos = nPos,
                                    r1NegFail = r1NegFail,r1PosFail = r1PosFail,
                                    rNegFail = rNegFail,rPosFail = rNegFail)
      
      np = ParallelTwoStageWorstFirst(pNeg = p0,pPos = p1, pr.pos = pr.pos,
                                    n1Neg = n1Neg,n1Pos = n1Pos,
                                    nNeg = nNeg,nPos = nPos,
                                    r1NegFail = r1NegFail,r1PosFail = r1PosFail,
                                    rNegFail = rNegFail,rPosFail = rNegFail)
      
      pp = ParallelTwoStageWorstFirst(pNeg = p1,pPos = p1, pr.pos = pr.pos,
                                    n1Neg = n1Neg,n1Pos = n1Pos,
                                    nNeg = nNeg,nPos = nPos,
                                    r1NegFail = r1NegFail,r1PosFail = r1PosFail,
                                    rNegFail = rNegFail,rPosFail = rNegFail)
      
      return(list(nn=nn,np=np,pp=pp))
    })
    
  } else if(approach == "bestFirst") {
    res = with(param,{
      nn = ParallelTwoStageBestFirst(pNeg = p0,pPos = p0, pr.pos = pr.pos,
                                      n1Neg = n1Neg,n1Pos = n1Pos,
                                      nNeg = nNeg,nPos = nPos,
                                      r1NegFail = r1NegFail,r1PosFail = r1PosFail,
                                      rNegFail = rNegFail,rPosFail = rNegFail)
      
      np = ParallelTwoStageBestFirst(pNeg = p0,pPos = p1, pr.pos = pr.pos,
                                      n1Neg = n1Neg,n1Pos = n1Pos,
                                      nNeg = nNeg,nPos = nPos,
                                      r1NegFail = r1NegFail,r1PosFail = r1PosFail,
                                      rNegFail = rNegFail,rPosFail = rNegFail)
      
      pp = ParallelTwoStageBestFirst(pNeg = p1,pPos = p1, pr.pos = pr.pos,
                                      n1Neg = n1Neg,n1Pos = n1Pos,
                                      nNeg = nNeg,nPos = nPos,
                                      r1NegFail = r1NegFail,r1PosFail = r1PosFail,
                                      rNegFail = rNegFail,rPosFail = rNegFail)
      
      return(list(nn=nn,np=np,pp=pp))
    })
    
  }
  
  dta = data.frame(matrix(0,0,13))
  names(dta) = c("Design","P(+)","P(b-)","P(b+)","interim stop","interim unselected","interim positive","Futility","Success (b- & b+)","Success (b+)","mean SS","99% SS","mean tested")
  
  nn = res$nn
  dta[1,] = c(approach,pr.pos,p0,p0,roundWZero(c(nn$intStop,nn$intBoth,nn$intPos),3),roundWZero(c(nn$h0,nn$h2,nn$h1),3),roundWZero(nn$ss,1),NA, roundWZero(nn$muScreen,1))
  
  np = res$np
  dta[2,] = c(approach,pr.pos,p0,p1,roundWZero(c(np$intStop,np$intBoth,np$intPos),3),roundWZero(c(np$h0,np$h2,np$h1),3),roundWZero(np$ss,1),NA, roundWZero(np$muScreen,1))
  
  pp = res$pp
  dta[3,] = c(approach,pr.pos,p1,p1,roundWZero(c(pp$intSto,pp$intBoth,pp$intPos),3),roundWZero(c(pp$h0,pp$h2,pp$h1),3),roundWZero(pp$ss,1),NA, roundWZero(pp$muScreen,1))
  
  return(dta)
}