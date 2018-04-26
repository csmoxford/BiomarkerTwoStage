


library(stats)
library(TailRank)

stopRuleEff = function(n,alpha,beta,threshold,probability) {
  res = which(sapply(1:n, function(r) 1-pbeta(threshold,r+alpha, n-r+beta)) > probability)
  if(length(res) == 0) {
    return(NA)
  }
  res[1]
}

stopRuleFut = function(n,alpha,beta,threshold,probability) {
  res = which(sapply(0:n, function(r) pbeta(threshold,r+alpha, n-r+beta)) > probability)
  if(length(res) == 0) {
    return(-1)
  }
  res[length(res)]-1
}


posteriorPredictivePower = function(rCur,nCur,alpha,beta,nEnd,sucessesNeeded) {
  if(rCur >= sucessesNeeded){
    return(1)
  }
  1-pbb(sucessesNeeded-rCur-1, nEnd-nCur, rCur+alpha,nCur-rCur+beta)
}



ParallelTwoStageBalancedStaggered = function(pNeg = 0.3, pPos = 0.3, pr.pos = 0.5, continueThresholdsPos, continueThresholdsNeg, continueThresholdsAll, n0Neg = 17, n0Pos = 17, nNeg = 34, nPos = 34) {

  # precompute probabilities
  negBinMatPos = matrix(1,35,35)
  negBinMatNeg = matrix(1,35,35)
  probMatPos = matrix(1,35,35)
  probMatNeg = matrix(1,35,35)
  for(i in 1:35){
    for(j in 1:35){
      probMatPos[i,j] = dbinom(j-1,i-1,pPos)
      probMatNeg[i,j] = dbinom(j-1,i-1,pNeg)
      negBinMatPos[i,j] = dnbinom(j-1,i,pr.pos)
      negBinMatNeg[i,j] = dnbinom(j-1,i,1-pr.pos)
    }
  }

  Both = 0
  Positive = 0
  Neither = 0
  ss = 0
  screen = 0 # extra patients beyond ss
  
  screenNeg = pr.pos/(1-pr.pos) # number of positives needed to see a negative (on average)
  screenPos = (1-pr.pos)/pr.pos # number of negatives needed to see a positive (on average)
  
  #############################################################
  # This is not an efficient implimentation but does the job
  #############################################################
  # negative first at both analyses
  total = 0
  for(x in 0:(n0Pos-1)) {
    for (y in 0:(nPos-1-x)){
      for(r0Neg in 0:n0Neg){
        for(r1Neg in 0:(nNeg-n0Neg)){
          for(r0Pos in 0:x) {
            for(r1Pos in 0:y) {
              for(r2Pos in 0:(nPos-x-y)) {
                
                prob = negBinMatNeg[n0Neg,x+1] * negBinMatNeg[nNeg-n0Neg,y+1] * probMatNeg[n0Neg+1,r0Neg+1] * probMatNeg[nNeg-n0Neg+1,r1Neg+1] * probMatPos[x+1,r0Pos+1] * probMatPos[y+1, r1Pos+1] * probMatPos[nPos-x-y+1,r2Pos+1]
                
                total = total + prob
                
                if(prob > 0) {
                if(r0Neg + r0Pos >= continueThresholdsAll[x+n0Neg+1] & r0Neg >= continueThresholdsNeg[n0Neg+1]) { # ok in all go to next analysis
                  
                  if(r0Neg + r0Pos + r1Neg + r1Pos >=continueThresholdsAll[x+y+nNeg+1] & r0Neg + r1Neg >= continueThresholdsNeg[nNeg+1]) { # ok in all stop trial for success
                    Both = Both + prob
                    ss = ss + prob*(x+y+nNeg)
                  } else if(r0Pos + r1Pos >= continueThresholdsPos[x+y+1]) { # ok in positive => go to final analysis
                    if(r0Pos + r1Pos + r2Pos >= continueThresholdsPos[nPos+1]) { # ok in positive at final analysis
                      Positive = Positive + prob
                      ss =ss + prob*(nNeg+nPos)
                      screen = screen + screenPos*prob*(nPos-x-y)
                    } else {
                      Neither = Neither + prob
                      ss =ss + prob*(nNeg+nPos)
                      screen = screen + screenPos*prob*(nPos-x-y)
                    }
                  } else { # stop at interim
                    Neither = Neither + prob
                    ss =ss + prob*(x+y+nNeg)
                    
                  }
                  
                  
                } else if(r0Pos >= continueThresholdsPos[x+1]) { # ok in positive => go to final analysis
                  if(r0Pos + r1Pos + r2Pos >= continueThresholdsPos[nPos+1]) { # ok in positive at final analysis
                    Positive = Positive + prob
                    ss =ss + prob*(nPos+n0Neg)
                    screen = screen + screenPos*prob*(nPos-x)
                  } else {
                    Neither = Neither + prob
                    ss =ss + prob*(nPos+n0Neg)
                    screen = screen + screenPos*prob*(nPos-x)
                  }
                  
                } else { # stop at interim
                  Neither = Neither + prob
                  ss =ss + prob*(n0Pos+n0Neg)
                }
                
                }
                
              }
            }
          }
        }
      }
    }
  }
    
  # positive first at both analyses
  for(x in 0:(n0Neg-1)) {
    for (y in 0:(nNeg-1-x)){
      for(r0Neg in 0:x){
        for(r1Neg in 0:y){
          for(r2Neg in 0:(nNeg-x-y)) {
            for(r0Pos in 0:n0Pos) {
              for(r1Pos in 0:(nPos-n0Pos)) {
              
                prob = negBinMatPos[n0Pos,x+1] * negBinMatPos[nPos-n0Pos,y+1] * probMatNeg[x+1, r0Neg+1] * probMatNeg[y+1, r1Neg+1] * probMatNeg[nNeg-x-y+1,r2Neg+1] * probMatPos[n0Pos+1,r0Pos+1] * probMatPos[nPos-n0Pos+1,r1Pos+1]
                
                total = total + prob
                
                if(prob > 0){
                
                if(r0Neg + r0Pos >= continueThresholdsAll[x+n0Pos+1] & r0Neg >= continueThresholdsNeg[x+1]) { # ok in all go to next analysis
                  
                  if(r0Neg + r0Pos + r1Neg + r1Pos >=continueThresholdsAll[x+y+nPos+1] & r0Neg + r1Neg >= continueThresholdsNeg[x+y+1]) { # ok in all continue to final analysis
                    if(r0Neg + r0Pos + r1Neg + r1Pos + r2Neg >= continueThresholdsAll[nNeg+nPos+1] & r0Neg + r1Neg + r2Neg >= continueThresholdsNeg[nNeg+1]) { # ok at final analysis
                      Both = Both + prob
                      ss = ss + prob*(nPos+nNeg)
                      screen = screen + screenNeg*prob*(nNeg-x-y)
                    } else if(r0Pos + r1Pos >= continueThresholdsPos[nPos+1]) { # positive only
                      Positive = Positive + prob
                      ss = ss + prob*(nPos+nNeg)
                      screen = screen + screenNeg*prob*(nNeg-x-y)
                    } else { # neither
                      Neither = Neither + prob
                      ss = ss + prob*(nPos+nNeg)
                      screen = screen + screenNeg*prob*(nNeg-x-y)
                    }
                    
                  } else if(r0Pos + r1Pos >= continueThresholdsPos[nPos+1]) { # ok in positive stop
                    Positive = Positive + prob
                    ss = ss + prob*(nPos+x+y)
                  } else { # stop at interim
                    Neither = Neither + prob
                    ss = ss + prob*(nPos+x+y)
                  }
                  
                  
                } else if(r0Pos >= continueThresholdsPos[n0Pos+1]) { # ok in positive => go to final analysis
                  
                  if(r0Pos + r1Pos >= continueThresholdsPos[nPos+1]) { # ok in positive at final analysis
                    Positive = Positive + prob
                    ss = ss + prob*(nPos+x)
                    screen = screen + screenPos*prob*(nPos-n0Pos)
                  } else {
                    Neither = Neither + prob
                    ss = ss + prob*(nPos+x)
                    screen = screen + screenPos*prob*(nPos-n0Pos)
                  }
                  
                } else { # stop at interim
                  Neither = Neither + prob
                  ss = ss + prob*(n0Pos+x)
                }
                }
                
              }
            }
          }
        }
      }
    }
  }
  
  
  # negative first at first analysis then positive
  for(x in 0:(n0Pos-1)) {
    for (y in 0:(nNeg-n0Neg-1)){
      for(r0Neg in 0:n0Neg){
        for(r1Neg in 0:(y)){
          for(r2Neg in 0:(nNeg-n0Neg-y)) {
            for(r0Pos in 0:x) {
              for(r1Pos in 0:(nPos-x)) {
                prob = negBinMatNeg[n0Neg,x+1] * negBinMatPos[nPos-x,y+1] * probMatNeg[n0Neg+1,r0Neg+1] * probMatNeg[y+1, r1Neg+1] * probMatNeg[nNeg-n0Neg-y+1,r2Neg+1] * probMatPos[x+1,r0Pos+1] * probMatPos[nPos - x+1,r1Pos+1]
                
                total = total + prob
                
                if(prob > 0) {
                if(r0Neg + r0Pos >= continueThresholdsAll[x+n0Neg+1] & r0Neg >= continueThresholdsNeg[n0Neg+1]) { # ok in all go to next analysis
                  
                  if(r0Neg + r0Pos + r1Neg + r1Pos >=continueThresholdsAll[n0Neg+y + nPos+1] & r0Neg + r1Neg >= continueThresholdsNeg[n0Neg+y+1]) { # ok in got to final analysis
                    if(r0Neg + r0Pos + r1Neg + r1Pos >=continueThresholdsAll[n0Neg+y + nPos+1] & r0Neg + r1Neg + r2Neg >= continueThresholdsNeg[nNeg+1]) {
                      Both = Both + prob
                      ss = ss + prob*(nPos+nNeg)
                      screen = screen + screenNeg*prob*(nNeg-n0Neg-y)
                    } else if(r0Pos + r1Pos >= continueThresholdsPos[x+y+1]) { # ok in positive stop
                      Positive = Positive + prob
                      ss = ss + prob*(nPos+nNeg)
                      screen = screen + screenNeg*prob*(nNeg-n0Neg-y)
                    } else { # bad in both
                      Neither = Neither + prob
                      ss = ss + prob*(nPos+nNeg)
                      screen = screen + screenNeg*prob*(nNeg-n0Neg-y)
                    }
                    
                  } else if(r0Pos + r1Pos >= continueThresholdsPos[nPos+1]) {
                    Positive = Positive + prob
                    ss = ss + prob*(nPos+n0Neg + y)
                  } else {
                    Neither = Neither + prob
                    ss = ss + prob*(nPos + n0Neg + y)
                  }
                  
                } else if(r0Pos >= continueThresholdsPos[x+1]) { # ok in positive => go to final analysis
                  
                  if(r0Pos + r1Pos >= continueThresholdsPos[nPos+1]) { # ok in positive at final analysis
                    Positive = Positive + prob
                    ss = ss + prob*(nPos + n0Neg)
                    screen = screen + screenPos*prob*(nPos-x)
                  } else {
                    Neither = Neither + prob
                    ss = ss + prob*(nPos + n0Neg)
                    screen = screen + screenPos*prob*(nPos-x)
                  }
                  
                } else { # stop at interim
                  Neither = Neither + prob
                  ss = ss + prob*(x + n0Neg)
                }
                }
              }
            }
          }
        }
      }
    }
  }
  
  
  # positive first at first analysis then negative
  for(x in 0:(n0Neg-1)) {
    for (y in 0:(nPos-n0Pos-1)){
      for(r0Neg in 0:x){
        for(r1Neg in 0:(nNeg-x)){
          for(r0Pos in 0:n0Pos) {
            for(r1Pos in 0:(y)) {
              for(r2Pos in 0:(nPos-n0Pos-y)) {
                
                prob = negBinMatPos[n0Pos,x+1] * negBinMatNeg[nNeg-x,y+1] * probMatNeg[x+1,r0Neg+1] * probMatNeg[nNeg - x +1, r1Neg+1] * probMatPos[n0Pos+1,r0Pos+1] * probMatPos[y +1,r1Pos+1] * probMatPos[nPos - n0Pos - y +1,r2Pos+1]
                
                total = total + prob
          
                
                if(prob > 0) {
                if(r0Neg + r0Pos >= continueThresholdsAll[x+n0Pos+1] & r0Neg >= continueThresholdsNeg[x+1]) { # ok in all go to next analysis
                  
                  if(r0Neg + r0Pos + r1Neg + r1Pos >=continueThresholdsAll[nNeg + n0Pos + y+1] & r0Neg + r1Neg >= continueThresholdsNeg[nNeg+1]) { # ok in all continue declare both positive
                    Both = Both + prob
                    ss = ss + prob*(n0Pos+y+nNeg)
                  } else if(r0Pos + r1Pos >= continueThresholdsPos[x+y+1]) {
                    if(r0Neg + r0Pos + r1Neg + r1Pos + r2Pos >= continueThresholdsAll[nNeg+nPos+1] & r0Neg + r1Neg >= continueThresholdsNeg[nNeg+1]) { # ok at final analysis
                      Both = Both + prob
                      ss = ss + prob*(nPos+nNeg)
                      screen = screen + screenPos*prob*(nPos-n0Pos-y)
                    } else if(r0Pos + r1Pos + r2Pos >= continueThresholdsPos[nPos+1]) { # positive only
                      Positive = Positive + prob
                      ss = ss + prob*(nPos+nNeg)
                      screen = screen + screenPos*prob*(nPos-n0Pos-y)
                    } else { # neither
                      Neither = Neither + prob
                      ss = ss + prob*(nPos+nNeg)
                      screen = screen + screenPos*prob*(nPos-n0Pos-y)
                    }
                    
                  }  else { # stop at interim
                    Neither = Neither + prob
                    ss = ss + prob*(n0Pos+y+nNeg)
                  }
                  
                  
                } else if(r0Pos >= continueThresholdsPos[n0Pos+1]) { # ok in positive => go to final analysis
                  
                  if(r0Pos + r1Pos + r2Pos >= continueThresholdsPos[nPos+1]) { # ok in positive at final analysis
                    Positive = Positive + prob
                    ss = ss + prob*(nPos+x)
                    screen = screen + screenPos*prob*(nPos-n0Pos)
                  } else {
                    Neither = Neither + prob
                    ss = ss + prob*(nPos+x)
                    screen = screen + screenPos*prob*(nPos-n0Pos)
                  }
                  
                } else { # stop at interim
                  Neither = Neither + prob
                  ss = ss + prob*(n0Pos+x)
                }
                
                }
              }
            }
          }
        }
      }
    }
  }
  
  screen = screen + ss
  print(total)
  return(list(Both=Both, Positive=Positive, Neither=Neither, ss = ss, muScreen = screen))
}


# continueThresholdsPos, continueThresholdsNeg, continueThresholdsAll are the number of responses required to continue if there are i-1 patients (zero indexed)
ParallelTwoStageBalancedStaggeredProperties = function(p0 = 0.1,p1 = 0.3, pr.pos = 0.5,
                                                       n1Neg = 17,n1Pos = 17,
                                                                  nNeg = 34,nPos = 34,
                                                                  continueThresholdsPos, continueThresholdsNeg, continueThresholdsAll) {
  
  
  
  nn = ParallelTwoStageBalancedStaggered(pNeg = p0, pPos = p0, pr.pos = pr.pos, continueThresholdsPos, continueThresholdsNeg, continueThresholdsAll, n0Neg = n1Neg, n0Pos = n1Pos, nNeg = nNeg, nPos = nPos)
  
  print(nn)
  
  np = ParallelTwoStageBalancedStaggered(pNeg = p0, pPos = p1, pr.pos = pr.pos, continueThresholdsPos, continueThresholdsNeg, continueThresholdsAll, n0Neg = n1Neg, n0Pos = n1Pos, nNeg = nNeg, nPos = nPos)
  
  print(np)
  
  pp = ParallelTwoStageBalancedStaggered(pNeg = p1, pPos = p1, pr.pos = pr.pos, continueThresholdsPos, continueThresholdsNeg, continueThresholdsAll, n0Neg = n1Neg, n0Pos = n1Pos, nNeg = nNeg, nPos = nPos)
  
  print(pp)
  
  
  dta = data.frame(matrix(0,0,13))
  names(dta) = c("Design","P(+)","P(b-)","P(b+)","interim stop","interim unselected","interim positive","Futility","Success (b- & b+)","Success (b+)","mean SS","99% SS","mean tested")
  
  dta[1,] = c("Balanced Staggered",pr.pos,p0,p0,NA,NA,NA,roundWZero(c(nn$Neither,nn$Positive,nn$Both),3),roundWZero(nn$ss,1),NA, roundWZero(nn$muScreen,1))
  
  dta[2,] = c("Balanced Staggered",pr.pos,p0,p1,NA,NA,NA,roundWZero(c(np$Neither,np$Positive,np$Both),3),roundWZero(np$ss,1),NA, roundWZero(np$muScreen,1))
  
  dta[3,] = c("Balanced Staggered",pr.pos,p1,p1,NA,NA,NA,roundWZero(c(pp$Neither,pp$Positive,pp$Both),3),roundWZero(pp$ss,1),NA, roundWZero(pp$muScreen,1))
  
  return(dta)
  
  
}
