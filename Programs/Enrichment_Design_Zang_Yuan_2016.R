#######################################################################################################################
##  Optimal sequential enrichment designs for phase II clinical trials
##          by Yong Zang and Ying Yuan
##
## This file contains the function for implementation of the proposed optimal and minimax sequential designs
##
########################################################################################################################

library(clinfun)

seq = function(p0,p1,alpha1,power1,alpha0,power0,n.max,u,method) {
  ## p0: unacceptable response rate
  ## p1: desirable response rate
  ## alpha1: type I error for marker-positive subgroup
  ## power1: power for marker-positive subgroup
  ## alpha0: type I error for marker-negative subgroup
  ## power0: power for marker-negative subgroup
  ## n.max: maximum sample size for each marker subgroup
  ## u: the upper bound of response rate for marker positive subgroup
  ## method: "opt" refers to the optimal design; "min" refers to the minmax design
  beta1=1-power1 ## type II error for  marker-positive group
  prob=function(r1,n1,r,n,ppu){
    a=rep(0, (n1-r1))
    for (i in 1 : n1-r1){
      a[i]=dbinom(i+r1,n1,ppu)*(1-pbinom(r-i-r1,n-n1,ppu))

    }
    return(sum(a))
  }
  if(method=="opt"){
    pos=ph2simon(p0,p1,alpha1,beta1,n.max)
    indexp=order(pos[[5]][,5])[1]
    desp=pos[[5]][indexp,][1:4] ## design parameters for marker-positive subgroup
    r1p=desp[[1]]
    n1p=desp[[2]]
    rp=desp[[3]]
    np=desp[[4]]
    p=prob(r1p,n1p,rp,np,u)
    alpha0=alpha0/p ## adjusted tyep I error for marker-negative group
    beta0=(power1-power0)/power1 ## adjusted type II error for marker-negative group
    neg=ph2simon(p0,p1,alpha0,beta0,n.max)
    indexn=order(neg[[5]][,5])[1]
    desn=neg[[5]][indexn,][1:4] ## design parameters for marker-negative subgroup
    r1n=desn[[1]]
    n1n=desn[[2]]
    rn=desn[[3]]
    nn=desn[[4]]
    pet1=pbinom(r1p,n1p,p0)
    pet2=0
    for (i in (r1p+1): n1p){
      a=dbinom(i, n1p, p0)*(1-pbinom((rp-i),(np-n1p),p0))
      pet2=a+pet2
    }
    pet2=1-pet2
    pet3=1-(1-pet2)*(1-pbinom(r1n,n1n,p0))
    en=n1p+(1-pet1)*(np-n1p)+(1-pet2)*n1n+(1-pet3)*(nn-n1n)
    return(list("r1.positive"=desp[[1]],"n1.positve"=desp[[2]],"r.positive"=desp[[3]],"n.positive"=desp[[4]],"r1.negative"=desn[[1]],"n1.negative"=desn[[2]],"r.negative"=desn[[3]],"n.negative"=desn[[4]],"EN(p0)"=en,"PET1(p0)"=pet1,"PET2(p0)"=pet2,"PET3(p0)"=pet3))
  }
  if(method=="min"){
    pos=ph2simon(p0,p1,alpha1,beta1,n.max)
    indexp=order(pos[[5]][,4])[1]
    desp=pos[[5]][indexp,][1:4] ## design parameters for marker-positive subgroup
    r1p=desp[[1]]
    n1p=desp[[2]]
    rp=desp[[3]]
    np=desp[[4]]
    p=prob(r1p,n1p,rp,np,u)
    alpha0=alpha0/p ## adjusted tyep I error for marker-negative group
    beta0=(power1-power0)/power1 ## adjusted type II error for marker-negative group
    neg=ph2simon(p0,p1,alpha0,beta0,n.max)
    indexn=order(neg[[5]][,4])[1]
    desn=neg[[5]][indexn,][1:4] ## design parameters for marker-negative subgroup
    r1n=desn[[1]]
    n1n=desn[[2]]
    rn=desn[[3]]
    nn=desn[[4]]
    pet1=pbinom(r1p,n1p,p0)
    pet2=0
    for (i in (r1p+1): n1p){
      a=dbinom(i, n1p, p0)*(1-pbinom((rp-i),(np-n1p),p0))
      pet2=a+pet2
    }
    pet2=1-pet2
    pet3=1-(1-pet2)*(1-pbinom(r1n,n1n,p0))
    en=n1p+(1-pet1)*(np-n1p)+(1-pet2)*n1n+(1-pet3)*(nn-n1n)
    return(list("r1.positive"=desp[[1]],"n1.positve"=desp[[2]],"r.positive"=desp[[3]],"n.positive"=desp[[4]],"r1.negative"=desn[[1]],"n1.negative"=desn[[2]],"r.negative"=desn[[3]],"n.negative"=desn[[4]],"EN(p0)"=en,"PET1(p0)"=pet1,"PET2(p0)"=pet2,"PET3(p0)"=pet3))
  }
}



