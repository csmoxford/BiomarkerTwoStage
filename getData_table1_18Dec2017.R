
##############################################################################
# Independent parallel
source(paste0(programDir, 'Programs/Exact/IndependentParallel_GetExactProperties.R'))

p0=0.1
p1=0.3
pr.positive = 0.5
n1 = 17
n = 34

eta=c(0.98,0.95)
zeta=c(NA,0.90)
prior.a = 0.2/4
prior.b = 0.8/4

df = getIndParallelProperties("IndParallel",pr.positive, p0,p1, res, n1, n, eta, zeta, prior.a,prior.b)



##############################################################################
# Tandem Two-Stage
source(paste0(programDir, 'Programs/Exact/TandemTwoStage_GetExactProperties.R'))

dta = tandemTwoStageproperties(p0,p1,pr.pos = 0.5)
df = rbind(df,dta)

##############################################################################
# Parallel Two-stage
source(paste0(programDir, 'Programs/Exact/ParallelTwoStage_GetExactProperties.R'))

## WorstFirst
dta = ParallelTwoStageProperties(approach = "worstFirst",p0,p1, pr.pos = 0.5)
df = rbind(df,dta)

## BetsFirst
dta = ParallelTwoStageProperties(approach = "bestFirst",p0,p1, pr.pos = 0.5)
df = rbind(df,dta)

## Balanced
dta = ParallelTwoStageProperties(approach = "balanced",p0,p1, pr.pos = 0.5)
df = rbind(df,dta)
##############################################################################
# Sequential enrichment

source(paste0(programDir, 'Programs/Exact/SequentialEnrichment_GetExactProperties.R'))
library(EurosarcBayes)

resPositive = bayes_binom_one_postprob_nstage(c(17,34),p0=0.1, p1=0.3, eta=c(0.98,0.95), zeta=c(NA,0.90), 0.2/4, 0.8/4, round=TRUE)

resNegative = bayes_binom_one_postprob_nstage(c(17,34),p0=0.1, p1=0.3, eta=c(0.98,0.95), zeta=c(NA,0.90), 0.2/4, 0.8/4, round=TRUE)

GetPropertiesBayesianEnrichment = function(resPositive, resNegative, designName =  "Sequential Enrichment", pr.positive = 0.5) {
  
  design = list(
    n1.positve = resPositive@reviews[1],
    r1.positive = resPositive@failure[1],
    n.positive = resPositive@reviews[2],
    r.positive = resPositive@failure[2],
    n1.negative = resNegative@reviews[1],
    r1.negative = resNegative@failure[1],
    n.negative = resNegative@reviews[2],
    r.negative = resNegative@failure[2]
  )
  
  return(enrichGetPropertiesTable(design,designName, pr.positive))
}
dta = GetPropertiesBayesianEnrichment(resPositive, resNegative)

df = rbind(df,dta)


load(paste0(programDir,"BalancedStaggeredData_02Feb2018.Rdata"))
df = rbind(df, balancedData[balancedData$`P(+)` == "0.5",])
