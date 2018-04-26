library(EurosarcBayes)

# Bayes sample size


bayes_binom_one_postprob_onestage(0.1, 0.3, eta = 0.95, zeta = 0.95, 0.2/4, 0.8/4, round=TRUE) # 34
bayes_binom_one_postprob_nstage(c(34),p0=0.1, p1=0.3, eta=c(0.95), zeta=c(0.95), 0.2/4, 0.8/4, round=TRUE) # 34

# with interim analysis
bayes_binom_one_postprob_nstage(c(17,34),p0=0.1, p1=0.3, eta=c(0.98,0.95), zeta=c(0.98,0.90), 0.2/4, 0.8/4, round=TRUE)
# failure criteria at interim is 1 or fewer
# failure criteria at final analysis is 6 or fewer
# By design we use the same properties for the positive and negative groups. This need not be the case.


##################################################
# in the balanced design we also check the unselected group prior to looking at the groups individually.
# fixed interim with 17 in each group
resCombined = bayes_binom_one_postprob_nstage(c(34),p0=0.1, p1=0.3, eta=c(0.98), zeta=c(0.95), 0.2/4, 0.8/4, round=TRUE)
# failure criteria is 5 or fewer

# fixed final analysis with 34 in each group
resCombined = bayes_binom_one_postprob_nstage(c(68),p0=0.1, p1=0.3, eta=c(0.95), zeta=c(0.95), 0.2/4, 0.8/4, round=TRUE)
# success criteria is 12 (failure is 11)
