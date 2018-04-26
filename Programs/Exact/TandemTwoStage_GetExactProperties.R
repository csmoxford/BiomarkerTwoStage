
tandemTwoStageproperties = function(
  p0=0.1,p1=0.3, pr.pos = 0.5,
  n1 = 17,n = 34,
  failInt = 1,failEnd = 6) {
  
  dta = data.frame(matrix(0,0,13))
  names(dta) = c("Design","P(+)","P(b-)","P(b+)","interim stop","interim unselected","interim positive","Futility","Success (b- & b+)","Success (b+)","mean SS","99% SS","mean tested")
  
  nn = tandemTwoStageExact(p0,p0,pr.pos,n1,n,failInt,failEnd)
  dta[1,] = c("tandemTwoStage",pr.pos,p0,p0,roundWZero(c(nn$intStop,nn$intBoth,nn$intPos),3),roundWZero(c(nn$h0,nn$h2,nn$h1),3),roundWZero(nn$ss,1),NA, roundWZero(nn$muScreen,1))
  
  nn = tandemTwoStageExact(p0,p1,pr.pos,n1,n,failInt,failEnd)
  dta[2,] = c("tandemTwoStage",pr.pos,p0,p1,roundWZero(c(nn$intStop,nn$intBoth,nn$intPos),3),roundWZero(c(nn$h0,nn$h2,nn$h1),3),roundWZero(nn$ss,1),NA, roundWZero(nn$muScreen,1))
  
  nn = tandemTwoStageExact(p1,p1,pr.pos,n1,n,failInt,failEnd)
  dta[3,] = c("tandemTwoStage",pr.pos,p1,p1,roundWZero(c(nn$intStop,nn$intBoth,nn$intPos),3),roundWZero(c(nn$h0,nn$h2,nn$h1),3),roundWZero(nn$ss,1),NA, roundWZero(nn$muScreen,1))
  
  return(dta)
}

tandemTwoStageExact = function(
  pNeg = 0.1,pPos = 0.3, pr.pos = 0.5,
  n1 = 17,n = 34,
  failInt = 1,failEnd = 6
) {
  
  
  prSuccessBoth = 0
  prSuccessPos = 0
  prFailBoth = 0
  prFailPos = 0
  prFailPosInt = 0
  # total = 0
  
  ssExtra = 0
  
  
  pnPos = dbinom(0:n1,n1,pr.pos)
  for(n1Pos in 0:17) {
    
    r1Pos = dbinom(0:n1Pos,n1Pos,pPos) 
    r1Neg = dbinom(0:(n1-n1Pos),(n1-n1Pos), pNeg)
    
    r1 = rep(0,18) # continue (distribution of response if continuing unselected)
    r1PosF = rep(0,length(r1Pos)) # continue +ve only (distribution of response in positive if +ve only)
    for(i in 1:length(r1Pos)) {
      for(j in 1:length(r1Neg)) {
        prob = r1Pos[i]*r1Neg[j]*pnPos[n1Pos+1]
        if(i+j-2 <= failInt) {
          r1PosF[i] = r1PosF[i] + prob
        } else {
          r1[i+j-1] = r1[i+j-1] + prob
        }
      }
    }
    # total = total + sum(r1PosF)+sum(r1)
    # continue both
    pnPos2 = dbinom(0:(n-n1),(n-n1),pr.pos)
    for(nPos in 0:17) {
      rPos = dbinom(0:nPos,nPos,pPos) 
      rNeg = dbinom(0:(n-n1-nPos),(n-n1-nPos), pNeg)
      for(i in 1:length(rPos)) {
        for(j in 1:length(rNeg)) {
          # print(paste0(i+j-2,", ",max(1,(failEnd+1-i-j+2))))
          prSuccessBoth = prSuccessBoth + rPos[i]*rNeg[j]*pnPos2[nPos+1]*sum(r1[max(1,(failEnd-i-j+4)):(n1+1)])
          if(i+j-1 <= failEnd) {
            prFailBoth = prFailBoth + rPos[i]*rNeg[j]*pnPos2[nPos+1]*sum(r1[1:(failEnd-i-j+3)])
          }
        }
      }
    }
    
    # continue +ve only
    # Get to interim first
    pr = dbinom(0:(n1-n1Pos),(n1-n1Pos),pPos)
    r1PosInt = rep(0,n1+1)
    for(i in 1:length(r1PosF)) {
      for(j in 1:length(pr)) {
        r1PosInt[i+j-1] = r1PosInt[i+j-1] + r1PosF[i]*pr[j]
        ssExtra = ssExtra + r1PosF[i]*pr[j]*(n1-n1Pos)
      }
    }
    
    prFailPosInt = prFailPosInt + sum(r1PosInt[1:2])
    r1PosInt[1:2] = 0
    
    rPos = dbinom(0:(n-n1),(n-n1),pPos)
    for(i in 1:length(r1PosInt)) {
      prSuccessPos = prSuccessPos + r1PosInt[i]*sum(rPos[max(1,(failEnd+1-i+2)):(n-n1+1)])
      if(i <= failEnd) {
        prFailPos = prFailPos + r1PosInt[i]*sum(rPos[1:(failEnd-i+2)])
      }
    }
    
  }
  
  h0 = prFailPosInt + prFailBoth + prFailPos
  h1 = prSuccessPos
  h2 = prSuccessBoth
  
  ss = 17 + ssExtra + 17*(prFailPos+prFailBoth+prSuccessBoth+prSuccessPos)
  muScreen = 17 + ssExtra/pr.pos + 17*(prSuccessPos+prFailPos)/pr.pos + 17*(prSuccessBoth+prFailBoth)
  
  return(list(h0=h0, h1=h1, h2=h2,ss=ss, muScreen = muScreen, intStop = NA, intPos = NA,intBoth = NA))
}

return(list())