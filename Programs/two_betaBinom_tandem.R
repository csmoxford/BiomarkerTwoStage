
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

  # compute the response outcome
  # print(c(decision@stage, decision@part))
  if(decision@stage == "1" | (decision@stage == "2" & data$biomarker[rw] > 0)){
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
      if(nPat >= p$nInt[1]) {
        decision@analyse = TRUE
      }
    } else {
      if(nPat >= p$nFinal[1]) {
        decision@analyse = TRUE
      }
    }

  } else if(decision@stage == "2") {
    if(decision@part == "a") {
      if(nPat1 >= p$nInt[2]) {
        decision@analyse = TRUE
      }
    } else {
      if(nPat1 >= p$nFinal[2]) {
        decision@analyse = TRUE
      }
    }

  }

  return(decision)
})


theModel = experimentalModel(fun = function(data, ...){
  p = list(...) # target, doseSupport, MTDq


  ########################################
  # to use when not using complete data
  curTime = max(data$arrivalTime)

  #######################################
  # prepare data to fit jags model
  # patients with outcome
  data = data[!is.na(data$outcome),]

  post0 = c(p$prior0[1] + sum(data$outcome==1 & data$biomarker == 0), p$prior0[2] + sum(data$outcome==0 & data$biomarker == 0))
  post1 = c(p$prior1[1] + sum(data$outcome==1 & data$biomarker == 1), p$prior1[2] + sum(data$outcome==0 & data$biomarker == 1))

  postAll = post0 + post1

  return(list(post0 = post0, post1 = post1, postAll = postAll))
})




theDecision = makeDecisions(fun = function(data, model, ...) {
  p = list(...) # nMax, nDoseMax, doseSupport, target, MTDq, cohortSize

  postProb = function(prior, threshold, side = 0){
    pr = pbeta(threshold, prior[1], prior[2])
    if(side){
      pr = 1 - pr
    }
    return(pr)
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
      if(postProb(model$postAll,p$biomarkerH1.0, 0) < p$postProbIF){
        stage = "1"
        part = "b"
        trail = c(trail, "1b")

      } else {
        stage = "2"
        part = "a"
        trail = c(trail, "2a")
      }
    } else if(p$decision@stage == "2") {
      if(postProb(model$post1,p$biomarkerH1.1, 0) < p$postProbIF){
        stage = "2"
        part = "b"
        trail = c(trail, "2b")
      } else {
        continue = FALSE
        stopReason = "2a futility"
      }
    } else {
      stop("what stage is this?")
    }



  } else { # final
    continue = FALSE

    if(p$decision@stage == "2") {

      if(postProb(model$post1,p$biomarkerH0.1, 1) > p$postProb) {
        stopReason = "2b success 1"
      } else {
        stopReason = "2b futility"
      }

    } else if(p$decision@stage == "1") {

      if(postProb(model$postAll,p$biomarkerH0.0, 1) > p$postProb) {
        stopReason = "1b success 01"
      } else {
        stopReason = "1b futility"
      }

    } else {
      stop("what stage is this?")
    }


  }

  if(nPat > 400){
    continue = FALSE
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
    biomarkerResp = c(0.1, 0.1), # probability of response for each level
    biomarkerH0.0 = 0.1, # posterior probability for efficacy is posterior greater
    biomarkerH0.1 = 0.1, # than this value
    biomarkerH1.0 = 0.3, # posterior probability for futility is posterior lesss
    biomarkerH1.1 = 0.3, # than this value
    prior0 = c(0.2,0.8)/4, # Beta prior for b^n
    prior1 = c(0.2,0.8)/4, # Beta prior for b^p
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
}

