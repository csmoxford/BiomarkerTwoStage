library(EurosarcBayes)

# Bayes sample size


bayes_binom_one_postprob_onestage(0.1, 0.3, eta = 0.95, zeta = 0.95, 0.2/4, 0.8/4, round=TRUE) # 34
bayes_binom_one_postprob_nstage(c(34),p0=0.1, p1=0.3, eta=c(0.95), zeta=c(0.95), 0.2/4, 0.8/4, round=TRUE) # 34

# with interim analysis
bayes_binom_one_postprob_nstage(c(17,34),p0=0.1, p1=0.3, eta=c(0.98,0.95), zeta=c(0.98,0.90), 0.2/4, 0.8/4, round=TRUE)

