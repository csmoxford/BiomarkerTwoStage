
source(paste0(programDir, 'Programs/utilFunctions/enrichGetPropertiesTable.R'))
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
