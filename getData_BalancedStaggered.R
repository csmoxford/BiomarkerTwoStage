library(EurosarcBayes)
library(peter.utilities)
library(prettyTables)
source('I:/Data/Peter_Dutton/Papers/wee1Paper/BiomarkerTwoStage/Programs/Exact/ParallelTwoStage_staggered_GetExactProperties.R')




n = 34

# prior
alpha = 0.2
beta = 0.8

r = stopRuleEff(n, alpha, beta,0.1,0.95)
stopRuleFut(n, alpha, beta,0.3,0.95)


threshold = rep(NA,n)
for(ncur in 1:n) {
  threshold[ncur] = which(sapply( 0:ncur, function(rcur) posteriorPredictivePower(rcur,ncur,alpha,beta,n,r)) > 0.05)[1]-1
}
threshold = c(0,threshold)



rall = stopRuleEff(2*n, alpha, beta,0.1,0.95)


thresholdComb = rep(NA,2*n)
for(ncur in 1:(2*n)) {
  thresholdComb[ncur] = which(sapply( 0:ncur, function(rcur) posteriorPredictivePower(rcur,ncur,alpha,beta,2*n,rall)) > 0.05)[1]-1
}
thresholdComb = c(0,thresholdComb)

##############################################################################
prevalences = c(0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8)


dfnames = c("Design","P(+)","P(b-)","P(b+)","interim stop","interim unselected","interim positive","Futility","Success (b- & b+)","Success (b+)","mean SS","99% SS","mean tested")
df = data.frame(matrix(0,nrow=0,ncol=length(dfnames)), stringsAsFactors = FALSE)
names(df) = dfnames

##############################################################################


for(pr.positive in prevalences) {

  dta = ParallelTwoStageBalancedStaggeredProperties(
    p0 = 0.1,
    p1 = 0.3,
    pr.pos = pr.positive,
    continueThresholdsPos = threshold,
    continueThresholdsNeg = threshold,
    continueThresholdsAll = thresholdComb,
    n1Neg = 17,
    n1Pos = 17,
    nNeg = 34,
    nPos = 34
  )
  
  df = rbind(df,dta)
}


save(df, file = paste0(outPath,"Figure5_BalancedStaggered.rData"))