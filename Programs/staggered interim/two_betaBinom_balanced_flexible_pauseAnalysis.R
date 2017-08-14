
library(trialSim)

baseline = simulateBaseline(fun = function(data = data.frame(), ...){
  p = list(...)

  if(dim(data)[1] == 0){
    cnames = c("id", "arrivalTime", "modelSeed", "biomarker","outcome","outcomeTime", "cohort","tested")
    data = data.frame(matrix(rep(0,0),ncol=length(cnames)))
    colnames(data) = cnames
  }



  modelSeed = as.integer(floor(runif(1,1,10^9)))
  id = max(0, data$id, na.rm = TRUE) + 1
  arrivalTime = max(0, data$arrivalTime, na.rm = TRUE) + rexp(1, p$arrivalRate)

  # remove untested individuals they are not needed??
  # print(sum(is.na(data$tested)))
  dm = dim(data)
  rng = runif(1)
  biomarker = ifelse(rng<p$biomarkerProb[1],0,1)

  data[dm[1]+1,1:4] = c(id, arrivalTime, modelSeed, biomarker)
  dm[1] = dm[1] + 1

  if(p$decision@part == "both"){
    data$tested[dm[1]] = FALSE
  }

  return(data)
})


outcome = simulateOutcome(fun = function(data, decision, ...){
  p = list(...) # alpha, beta, toxicityMaxTime


  # row to populate
  rw = dim(data)[1]
  # if trial is recruiting assign patient to cohort and dose

  # compute the response outcome
  # print(c(rw, decision@stage, decision@part, decision@recruit))

  nPat0 = sum(!is.na(data$outcome) & data$biomarker == 0)
  nPat1 = sum(!is.na(data$outcome) & data$biomarker == 1)

  if(decision@recruit == "both"){
    data$tested[rw] = TRUE
    if((data$biomarker[rw] == 0 & nPat0 < p$nFinal[1]) | (data$biomarker[rw] == 1 & nPat1 < p$nFinal[2])){
      data$outcome[rw] = rbinom(1,1,p$biomarkerResp[data$biomarker[rw]+1])
      data$cohort[rw] = decision@stage
      data$outcomeTime[rw] = data$arrivalTime[rw] + p$eventTime
    }
  } else if(decision@recruit == "+ve"){
    data$tested[rw] = TRUE
    if(data$biomarker[rw] == 1 & nPat1 < p$nFinal[2]){
      data$outcome[rw] = rbinom(1,1,p$biomarkerResp[data$biomarker[rw]+1])
      data$cohort[rw] = decision@stage
      data$outcomeTime[rw] = data$arrivalTime[rw] + p$eventTime
    }
  } else if(decision@recruit == "-ve"){
    data$tested[rw] = TRUE
    if(data$biomarker[rw] == 0 & nPat0 < p$nFinal[1]){
      data$outcome[rw] = rbinom(1,1,p$biomarkerResp[data$biomarker[rw]+1])
      data$cohort[rw] = decision@stage
      data$outcomeTime[rw] = data$arrivalTime[rw] + p$eventTime
    }
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

  curTime = max(data$arrivalTime)

  nPatOutcome = sum(!is.na(data$outcome) & data$outcomeTime < curTime)
  nPatOutcome0 = sum(!is.na(data$outcome) & data$biomarker == 0  & data$outcomeTime < curTime)
  nPatOutcome1 = sum(!is.na(data$outcome) & data$biomarker == 1  & data$outcomeTime < curTime)

  # print(nPat0,nPat1,nPatOutcome0,nPatOutcome1)

  if(nPat1 >= p$nFinal[2]){
    if(decision@recruit == "both"){
      decision@recruit = "-ve"
    } else if(decision@recruit == "+ve"){
      decision@recruit = "none"
    }
  }

  if(nPat0 >= p$nFinal[1]){
    if(decision@recruit == "both"){
      decision@recruit = "+ve"
    } else if(decision@recruit == "-ve"){
      decision@recruit = "none"
    }
  }

  decision@analyse = FALSE
  if(decision@stage == "aa") {

    if(nPat0 >= p$nInt[1] & nPat1 >= p$nInt[2]) {
      decision@analyse = nPatOutcome0 >= p$nInt[1] & nPatOutcome1 >= p$nInt[2]
      decision@recruit = "none"

    } else if(nPat0 >= p$nInt[1]) {
      decision@recruit = "+ve"
      if(nPat1 + p$nWait < p$nInt[2]) {
        decision@analyse = nPatOutcome0 >= p$nInt[1]
        decision@recruit = "none"
      }
    } else if(nPat1 >= p$nInt[2]) {
      decision@recruit = "-ve"
      if(nPat0 + p$nWait < p$nInt[1]) {
        decision@analyse = nPatOutcome1 >= p$nInt[2]
        decision@recruit = "none"
      }
    }
  } else if(decision@stage == "ab"){

    if(nPat0 >= p$nInt[1] & nPat1 >= p$nFinal[2]) {
      decision@analyse = nPatOutcome0 >= p$nInt[1] & nPatOutcome1 >= p$nFinal[2]
      decision@recruit = "none"
    } else if(nPat0 >= p$nInt[1]) {
      decision@recruit = "+ve"
      if(nPat1 + p$nWait < p$nFinal[2]) {
        decision@analyse = nPatOutcome0 >= p$nInt[1]
        decision@recruit = "none"
      }
    } else if(nPat1 >= p$nFinal[2]) {
      decision@recruit = "-ve"
      if(nPat0 + p$nWait < p$nInt[1]) {
        decision@analyse = nPatOutcome1 >= p$nFinal[2]
        decision@recruit = "none"
      }
    }

  } else if(decision@stage == "ba"){

    if(nPat0 >= p$nFinal[1] & nPat1 >= p$nInt[2]) {
      decision@analyse = nPatOutcome0 >= p$nFinal[1] & nPatOutcome1 >= p$nInt[2]
      decision@recruit = "none"
    } else if(nPat0 >= p$nFinal[1]) {
      decision@recruit = "+ve"
      if(nPat1 + p$nWait < p$nInt[2]) {
        decision@analyse = nPatOutcome0 >= p$nFinal[1]
        decision@recruit = "none"
      }
    } else if(nPat1 >= p$nInt[2]) {
      decision@recruit = "-ve"
      if(nPat0 + p$nWait < p$nFinal[1]) {
        decision@analyse = nPatOutcome1 >= p$nInt[2]
        decision@recruit = "none"
      }
    }
  } else if(decision@stage == "bb"){

    if(nPat0 >= p$nFinal[1] & nPat1 >= p$nFinal[2]) {
      decision@analyse = nPatOutcome0 >= p$nFinal[1] & nPatOutcome1 >= p$nFinal[2]
      decision@recruit = "none"
    } else if(nPat0 >= p$nFinal[1]) {
      decision@recruit = "+ve"
      if(nPat1 + p$nWait < p$nFinal[2]) {
        decision@analyse = nPatOutcome0 >= p$nFinal[1]
        decision@recruit = "none"
      }
    } else if(nPat1 >= p$nFinal[2]) {
      decision@recruit = "-ve"
      if(nPat0 + p$nWait < p$nFinal[1]) {
        decision@analyse = nPatOutcome1 >= p$nFinal[2]
        decision@recruit = "none"
      }
    }

  } else if(decision@stage == ".a"){
    if(nPat1 >= p$nInt[2]) {
      decision@analyse = nPatOutcome1 >= p$nInt[2]
      decision@recruit = "none"
    }
  } else if(decision@stage == "a."){
    if(nPat0 >= p$nInt[1]) {
      decision@analyse = nPatOutcome0 >= p$nInt[1]
      decision@recruit = "none"
    }
  } else if(decision@stage == ".b"){
    if(nPat1 >= p$nFinal[2]) {
      decision@analyse = nPatOutcome1 >= p$nFinal[2]
      decision@recruit = "none"
    }
  } else if(decision@stage == "b."){
    if(nPat0 >= p$nFinal[1]) {
      decision@analyse = nPatOutcome0 >= p$nFinal[1]
      decision@recruit = "none"
    }
  } else {
    stop("no stage found")
  }
  if(nPat > 400){
    decision@analyse = TRUE
  }

  decision@analysisTime = curTime


  return(decision)
})


theModel = experimentalModel(fun = function(data, ...){
  p = list(...) # target, doseSupport, MTDq

  data = data[!is.na(data$outcome),]
  data = data[data$outcomeTime <= p$decision@analysisTime,]


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
  stage = ""
  part = ""
  trail = p$decision@trail
  stopReason = ""
  recruit = "both"

  nPat = sum(!is.na(data$outcome) & data$outcomeTime <= p$decision@analysisTime)
  nPat0 = sum(!is.na(data$outcome) & data$biomarker == 0 & data$outcomeTime <= p$decision@analysisTime)
  nPat1 = sum(!is.na(data$outcome) & data$biomarker == 1 & data$outcomeTime <= p$decision@analysisTime)
  # print(p$decision@stage)
  # print(p$decision@part)

  efficacyThreshold = p$postProb
  futilityThreshold = p$postProbIF

  testAll = (postProb(model$postAll,p$biomarkerH1.0, 0) < futilityThreshold) | (postProb(model$postAll,p$biomarkerH0.0, 1) > efficacyThreshold)
  testNegative = (nPat0 < p$nFinal[1] & postProb(model$post0,p$biomarkerH1.0, 0) < futilityThreshold) | (nPat0 >= p$nFinal[1] & postProb(model$post0,p$biomarkerH0.0, 1) > efficacyThreshold)
  testPositive = (nPat1 < p$nFinal[2] & postProb(model$post1,p$biomarkerH1.0, 0) < futilityThreshold) | (nPat1 >= p$nFinal[2] & postProb(model$post1,p$biomarkerH0.0, 1) > efficacyThreshold)



  staging0 = (nPat0 >= p$nFinal[1]) + (nPat0 >= p$nInt[1]) + 1
  staging0 = c("a","b",".")[staging0]
  staging1 = (nPat1 >= p$nFinal[2]) + (nPat1 >= p$nInt[2]) + 1
  staging1 = c("a","b",".")[staging1]

  if(testAll) {
    # No evidence of futility in unselected group

    if(nPat0 >= p$nInt[1]) {
      # Enough marker negative patients to test marker-negative group

      if(testNegative) {
        # No evidence of futility in negative group
        recruit = "both"
      } else {
        staging0 = "."
        recruit = "+ve"
      }
    } else {
      # Not Enought patients to test marker negative group
      recruit = "both"
    }

  } else {
    # Evidence of futility in unselected group


    if(nPat1 >= p$nInt[2]) {
      # Enough marker positive patients to test marker-positive group

      if(testPositive) {
        # No evidence of futility in positive group
        recruit = "+ve"
        staging0 = "."
      } else {
        # Evidence of futility in positive group
        recruit = "stop"
      }

    } else {
      # Not Enought patients to test marker-positive group
      recruit = "+ve"
      staging0 = "."
    }
  }



  if(recruit == "stop"){
    continue = FALSE
    stopReason = "futility"
  } else if(nPat0 >= p$nFinal[1] & nPat1 >= p$nFinal[2]) {
    continue = FALSE
    if(testAll) {
      if(testNegative){
        stopReason = "success 01"
      } else {
        stopReason = "success 1"
      }
    } else {
      if(testPositive){
        stopReason = "success 1"
      } else {
        stopReason = "futility"
      }
    }
  } else if(nPat0 >= p$nFinal[1]) {
    if(testAll) {
      if(testNegative){
        stopReason = "success 01"
        continue = FALSE
      }
    }
  } else if(nPat1 >= p$nFinal[2]) {
    if(testPositive){
      if(p$decision@stage %in% c(".a",".b","..") | staging0 == "."){
        stopReason = "success 1"
        continue = FALSE
      }
    } else {
      stopReason = "futility"
      continue = FALSE
    }
  }

  # prevent the trial from relapsing from dropping -ve group
  if(p$decision@stage %in% c(".a",".b","..")) {
    staging0 = "."
  }
  if(p$decision@stage %in% c("a.","b.","..")) {
    staging1 = "."
  }

  stage = paste0(staging0, staging1)

  if(nPat > 400){
    continue = FALSE
  }

  # cat("\n",continue, stopReason, stage)
  print(c(model$postAll, model$post0, model$post1))
  print(c(continue, nPat0, nPat1,stage,stopReason, testAll,testNegative,testPositive, recruit))

  return(stageDecisionTwo(continue = continue, recruiting = TRUE, recruit = recruit, stopReason = stopReason, stage = stage, part = part, trail = trail))

})




setup = splitDataDesign(
  simBaseline = baseline,
  simOutcome = outcome,
  triggerAnalysis = triggerAnalysis,
  model = theModel,
  nSim = 1L, # number of simulations
  decision = theDecision, # This and the above plug in the function to run this trial for trialSim
  p = list( # options to adjust trial
    arrivalRate = 0.25, # arrival rate for timings
    eventTime = 60, # event time + time to analyse
    biomarkerProb = c(0.5, 0.5), # probability of each biomarker level
    biomarkerResp = c(0.1, 0.3), # probability of response for each level
    biomarkerH0.0 = 0.1, # posterior probability for efficacy is posterior greater
    biomarkerH0.1 = 0.1, # than this value
    biomarkerH1.0 = 0.3, # posterior probability for futility is posterior lesss
    biomarkerH1.1 = 0.3, # than this value
    prior0 = c(0.2,0.8)/4, # Beta prior for b^n
    prior1 = c(0.2,0.8)/4, # Beta prior for b^p
    postProb = 0.95, # posterior probability required for success at end
    postProbIF = 0.98, # posterior probability for futility to drop biomarker level
    postProbIAll = 0.98, # posterior probability for futility at early look to drop biomarker level
    nInt = c(17,17), # b^n, b^p counts to do interim
    nFinal = c(34,34), # b^n, b^p counts to do final (unselected stage 2)
    nWait = 0,
    decision = stageDecisionTwo(stage = "aa", part = "", recruit = "both", analyse = FALSE, trail = "1aa"), # this says we want to start a stage 1a. You shouldn't need to change this.
    sourceSetupCode = paste0(programDir,'Programs/staggered interim/simulationSetup.R')
  )
)


setup = splitDataDesign(
  simBaseline = baseline,
  simOutcome = outcome,
  triggerAnalysis = triggerAnalysis,
  model = theModel,
  nSim = 1L, # number of simulations
  decision = theDecision, # This and the above plug in the function to run this trial for trialSim
  p = list( # options to adjust trial
    arrivalRate = 0.25, # arrival rate for timings
    eventTime = 60, # event time + time to analyse
    biomarkerProb = c(0.5, 0.5), # probability of each biomarker level
    biomarkerResp = c(0.1, 0.3), # probability of response for each level
    biomarkerH0.0 = 0.1, # posterior probability for efficacy is posterior greater
    biomarkerH0.1 = 0.1, # than this value
    biomarkerH1.0 = 0.3, # posterior probability for futility is posterior lesss
    biomarkerH1.1 = 0.3, # than this value
    prior0 = c(0.2,0.8)/4, # Beta prior for b^n
    prior1 = c(0.2,0.8)/4, # Beta prior for b^p
    postProb = 0.95, # posterior probability required for success at end
    postProbIF = 0.98, # posterior probability for futility to drop biomarker level
    postProbIAll = 0.98, # posterior probability for futility at early look to drop biomarker level
    nInt = c(17,17), # b^n, b^p counts to do interim
    nFinal = c(34,34), # b^n, b^p counts to do final (unselected stage 2)
    nWait = 0,
    decision = stageDecisionTwo(stage = "aa", part = "", recruit = "both", analyse = FALSE, trail = "1aa"), # this says we want to start a stage 1a. You shouldn't need to change this.
    sourceSetupCode = paste0(programDir,'Programs/staggered interim/simulationSetup.R')
  )
)

# for testing only
run = FALSE
if(run){
  library(knitr)

  setup@seed = integer(0)
  setup@p$nWait = 3
  setup@p$biomarkerProb = c(0.3,0.7)
  res = simTrial(setup, nSim = 100L, parallel = FALSE)
  getSummary(res)

  nrecruited = sapply(res@sims, function(x) return(sum(!is.na(x$data$outcome))))
  max(nrecruited)
  stopreason = sapply(1:1000,function(i) res@sims[[i]]$decision@stopReason)
  which(stopreason == "1b success 01")
  which(nrecruited > 68)

  table( res@sims[[1]]$data$biomarker,res@sims[[1]]$data$cohort,res@sims[[1]]$data$tested, useNA = "ifany")

  i = 2
  res@sims[[i]]$decision@stopReason
  table(res@sims[[i]]$data$biomarker,res@sims[[i]]$data$outcome)

  setup@data = res@sims[[i]]$data

  model = getModel(setup)

  postProb(model$post0,0.1,1)
}







