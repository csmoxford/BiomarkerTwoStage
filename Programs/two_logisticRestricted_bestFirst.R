
library(trialSim)

baseline = simulateBaseline(fun = function(data = data.frame(), ...){
  p = list(...)


  if(dim(data)[1] == 0){
    cnames = c("id", "arrivalTime", "modelSeed", "biomarker","outcome","outcomeTime", "cohort")
    data = data.frame(matrix(rep(0,0),ncol=length(cnames)))
    colnames(data) = cnames
  }

  dm = dim(data)

  modelSeed = as.integer(floor(runif(1,1,10^9)))
  id = max(0, data$id, na.rm = TRUE) + 1
  arrivalTime = max(0, data$arrivalTime, na.rm = TRUE) + rexp(1, p$arrivalRate)
  rng = runif(1)
  biomarker = ifelse(rng<p$biomarkerProb[1],0,1)

  data[dm[1]+1,1:4] = c(id, arrivalTime, modelSeed, biomarker)
  dm[1] = dm[1] + 1
  return(data)
})


outcome = simulateOutcome(fun = function(data, decision, ...){
  p = list(...) # alpha, beta, toxicityMaxTime


  # row to populate
  rw = dim(data)[1]
  # if trial is recruiting assign patient to cohort and dose
  nPat0 = sum(!is.na(data$outcome) & data$biomarker == 0)
  nPat1 = sum(!is.na(data$outcome) & data$biomarker == 1)
  
  # print(c(nPat0,nPat1))
  
  # compute the response outcome
  # print(c(decision@stage, decision@part))
  if((decision@stage == "1" & (decision@part == "a" & ((data$biomarker[rw] == 0 & nPat0 < p$nInt[1]) | (data$biomarker[rw] == 1 & nPat1 < p$nInt[2])) | 
                               decision@part == "b")  & ((data$biomarker[rw] == 0 & nPat0 < p$nFinal[1]) | (data$biomarker[rw] == 1 & nPat1 < p$nFinal[2]))) | 
     (decision@stage == "2" & data$biomarker[rw] > 0)){
    data$outcome[rw] = rbinom(1,1,p$biomarkerResp[data$biomarker[rw]+1])
    data$cohort[rw] = paste0(decision@stage,decision@part)
    data$outcomeTime[rw] = data$arrivalTime[rw] + 1
  }

  return(data)
})



triggerAnalysis = readyForAnalysis(fun = function(data, decision, ...){
  p = list(...)
  # decision@recruiting
  # print(decision)
  nPat = sum(!is.na(data$outcome))
  nPat0 = sum(!is.na(data$outcome) & data$biomarker == 0)
  nPat1 = sum(!is.na(data$outcome) & data$biomarker == 1)

  decision@analyse = FALSE
  if(decision@stage == "1") {

    if(decision@part == "a") {
      if(nPat0 >= p$nInt[1] & nPat1 >= p$nInt[2]) {
        decision@analyse = TRUE
      }
    } else {
      if(nPat0 >= p$nFinal[1] & nPat1 >= p$nFinal[2]) {
        decision@analyse = TRUE
      }
    }

  } else if(decision@stage == "2") {

    if(nPat1 >= p$nFinal[2]) {
      decision@analyse = TRUE
    }

  }

  return(decision)
})


theModel = experimentalModel(fun = function(data, ...){
  p = list(...) # target, doseSupport, MTDq
  library(R2jags)
  data = data[!is.na(data$outcome),]
print(data)
  y = data$outcome
  x = as.numeric(data$biomarker>0)
  n = length(y)
  a = p$priora
  b = p$priorb

  # assume full information
  jags.model=function() {
    for(i in 1:n) {
      y[i] ~ dbern(theta[i])

      logit(theta[i]) <- alpha + beta * x[i]
    }

    alpha ~ dnorm(a[1], a[2])
    beta ~ dgamma(b[1], b[2])
  }

  # hard code the MCMC implimentation and parameters this time
  jagsVar = jagsVariables(model.file=jags.model,nSample=5000,burnin=500,nThin=1,n.chains=1,parameters.to.save=c("alpha","beta"))


  model = jags(data=c("y","x","n","a","b"),
               n.chains=jagsVar@n.chains,
               model.file=jagsVar@model.file,
               parameters.to.save=jagsVar@parameters.to.save,
               n.burnin=jagsVar@burnin,
               n.iter=jagsVar@nSample+jagsVar@burnin,
               n.thin=jagsVar@nThin,
               DIC=FALSE)

  return(model)
})




theDecision = makeDecisions(fun = function(data, model, ...) {
  p = list(...) # nMax, nDoseMax, doseSupport, target, MTDq, cohortSize


  inv_logit = function(alpha) 1/(1+exp(-alpha))

  postProbDist = function(model, level){
    m = model$BUGSoutput$sims.list
    if(level == 0){
      return(inv_logit(m$alpha))
    } else if(level == 1){
      return(inv_logit(m$alpha + m$beta))
    } else {
      stop("What level?")
    }
  }

  postProb = function(model, level, threshold, side = 0) {
    pr = mean(postProbDist(model, level) < threshold)
    if(side == 0){
      return(pr)
    } else {
      return(1-pr)
    }
  }

  continue = TRUE
  stage = p$decision@stage
  part = p$decision@part
  trail = p$decision@trail
  stopReason = ""

  nPat = sum(!is.na(data$outcome))
  nPat1 = sum(!is.na(data$outcome) & data$biomarker == 1)
  nPat2 = sum(!is.na(data$outcome) & data$biomarker == 2)

  if(p$decision@part == "a") { # interim

    if(p$decision@stage == "1") {
      # futility for group 0 not met
      if((nPat1 < p$nFinal[2] & postProb(model, 1,p$biomarkerH1.1, 0) < p$postProbIF) | (nPat1 >= p$nFinal[2] & postProb(model, 1,p$biomarkerH0.1, 1) > p$postProb)){
        if(postProb(model, 0,p$biomarkerH1.0, 0) < p$postProbIF){
              stage = "1"
              part = "b"
              trail = c(trail, "1b")
            # futility for group 1 and (group 1 and 2) not met
        } else {
          if((nPat1 >= p$nFinal[2] & postProb(model, 1,p$biomarkerH0.1, 1) > p$postProb)){
            stopReason = "1a success 1"
            continue = FALSE
          }
            stage = "2"
            part = "b"
            trail = c(trail, "2")
        }
      } else {
        stopReason = "1a futility"
        continue = FALSE
      }
    } else {
      stop("what stage is this?")
    }



  } else { # final
    continue = FALSE

    if(p$decision@stage == "2") {

      if(postProb(model, 1, p$biomarkerH0.1, 1) > p$postProb) {
          stopReason = "2 success 1"
      } else {
        stopReason = "2 futility"
      }

    } else if(p$decision@stage == "1") {

      if(postProb(model, 1, p$biomarkerH0.1, 1) > p$postProb) {
        if(postProb(model, 0, p$biomarkerH0.0, 1) > p$postProb) {
          stopReason = "1b success 01"
        } else {
          stopReason = "1b success 1"
        }

      } else {
        stopReason = "1b futility"
      }

    } else {
      stop("what stage is this?")
    }
  }

  return(stageDecision(continue = continue, recruiting = TRUE, stopReason = stopReason, stage = stage, part = part, trail = trail))

})


setup = splitDataDesign(
  simBaseline = baseline,
  simOutcome = outcome,
  triggerAnalysis = triggerAnalysis,
  model = theModel,
  nSim = 1L, # number of simulations
  decision = theDecision, # This and the above plug in the function to run this trial for trialSim
  p = list( # options to adjust trial
    arrivalRate = 1, # arrival rate for timings
    biomarkerProb = c(0.5, 0.5), # probability of each biomarker level
    biomarkerResp = c(0.1, 0.3), # probability of response for each level
    biomarkerH0.0 = 0.1, # posterior probability for efficacy is posterior greater
    biomarkerH0.1 = 0.1, # than this value
    biomarkerH1.0 = 0.3, # posterior probability for futility is posterior lesss
    biomarkerH1.1 = 0.3, # than this value
    priora = c(-1.5,0.1), # Beta prior for b^n
    priorb = c(1.5,1.5), # Beta prior for b^p
    postProb = 0.95, # posterior probability required for success at end
    postProbIF = 0.98, # posterior probability for futility to drop biomarker level
    nInt = c(17,17), # b^n, b^p counts to do interim
    nFinal = c(34,34), # b^n, b^p counts to do final (unselected stage 2)
    decision = stageDecision(stage = "1", part = "a", analyse = FALSE, trail = "1a") # this says we want to start a stage 1a. You shouldn't need to change this.
  )
)

run = FALSE
if(run){
  library(knitr)
  res = simTrial(setup,nSim = 1000L, parallel = TRUE)
  getSummary(res)
setup@p$biomarkerResp = c(0.1,0.1)
setup@p$biomarkerResp = c(0.1,0.3)
setup@p$biomarkerResp = c(0.3,0.3)
res = simTrial(setup,nSim = 10L, parallel = TRUE)
  i=4
  table(res@sims[[i]]$data$biomarker,res@sims[[i]]$data$outcome == 1)
  setup@data = res@sims[[i]]$data
  model = getModel(setup)
  getDecision(setup,model)

  setup@p$priorb = c(1,4)
  setup@p$priora = c(-1.5,2)

  n= 10000
  hist(rgamma(n,setup@p$priorb[1],setup@p$priorb[2]))
  hist(inv_logit(rnorm(n,setup@p$priora[1],sqrt(1/setup@p$priora[2]))))
  hist(inv_logit(rnorm(n,setup@p$priora[1],sqrt(1/setup@p$priora[2])) + rgamma(n,setup@p$priorb[1],setup@p$priorb[2])))


  par(mfrow=c(2,2))
  hist(model$BUGSoutput$sims.list$alpha)
  hist(model$BUGSoutput$sims.list$beta)
  hist(inv_logit(model$BUGSoutput$sims.list$alpha))
  hist(inv_logit(model$BUGSoutput$sims.list$alpha+model$BUGSoutput$sims.list$beta))

  stopreason = sapply(1:100,function(i) res@sims[[i]]$decision@stopReason)

  interim = sapply(1:100,function(i) res@sims[[i]]$data@stopReason)

  which(stopreason == "1b success 01")

}
