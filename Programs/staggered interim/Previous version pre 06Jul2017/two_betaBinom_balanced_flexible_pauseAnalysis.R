
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
      decision@recruit = "+ve"
    } else if(decision@recruit == "+ve"){
      decision@recruit = "none"
    }
  }

  if(nPat0 >= p$nFinal[1]){
    if(decision@recruit == "both"){
      decision@recruit = "-ve"
    } else if(decision@recruit == "-ve"){
      decision@recruit = "none"
    }
  }

  decision@analyse = FALSE
  if(decision@stage == "1aa") {

    if(nPatOutcome0 >= p$nInt[1] & nPatOutcome1 >= p$nInt[2]) {
      decision@analyse = TRUE
      decision@part = "both"
    } else if(nPatOutcome0 >= p$nInt[1] & nPatOutcome1 + p$nWait < p$nInt[2]) {
        decision@analyse = nPatOutcome0 >= p$nInt[1]
        decision@part = "-ve"
    } else if(nPatOutcome1 >= p$nInt[2] & nPatOutcome0 + p$nWait < p$nInt[1]){
        decision@analyse = TRUE
        decision@part = "+ve"
    }
  } else if(decision@stage == "2ab"){

    if(nPatOutcome0 >= p$nInt[1] & nPatOutcome1 >= p$nFinal[2]) {
      decision@analyse = TRUE
      decision@part = "both"
    } else if(nPatOutcome0 >= p$nInt[1]  & nPatOutcome1 + p$nWait < p$nFinal[2]) {
        decision@analyse = TRUE
        decision@part = "-ve"
    } else if(nPatOutcome1 >= p$nFinal[2]  & nPatOutcome0 + p$nWait < p$nInt[1]) {
        decision@analyse = TRUE
        decision@part = "+ve"
    }
  } else if(decision@stage == "2ba"){

    if(nPatOutcome0 >= p$nFinal[1] & nPatOutcome1 >= p$nInt[2]) {
      decision@analyse = TRUE
      decision@part = "both"
    } else if(nPatOutcome0 >= p$nFinal[1]  & nPatOutcome1 + p$nWait < p$nInt[2]) {
        decision@analyse = TRUE
        decision@part = "-ve"
    } else if(nPatOutcome1 >= p$nInt[2]  & nPatOutcome0 + p$nWait < p$nFinal[1]) {
        decision@analyse = TRUE
        decision@part = "+ve"
    }
  } else if(decision@stage == "3bb"){
    if(nPatOutcome0 >= p$nFinal[1] & nPatOutcome1 >= p$nFinal[2]) {
      decision@analyse = TRUE
      decision@part = "both"
    } else if(nPatOutcome0 >= p$nFinal[1]  & nPatOutcome1 + p$nWait < p$nFinal[2]) {
        decision@analyse = TRUE
        decision@part = "-ve"
    } else if(nPatOutcome1 >= p$nFinal[2]  & nPatOutcome0 + p$nWait < p$nFinal[1]) {
        decision@analyse = TRUE
        decision@part = "+ve"
    }
  } else if(decision@stage == "3fa"){
    if(nPatOutcome1 >= p$nInt[2]) {
      decision@analyse = TRUE
      decision@part = "+ve"
    }
  } else if(decision@stage == "3as"){
    if(nPatOutcome0 >= p$nInt[1]) {
      decision@analyse = TRUE
      decision@part = "-ve"
    }
  } else if(decision@stage == "4fb"){
    if(nPatOutcome1 >= p$nFinal[2]) {
      decision@analyse = TRUE
      decision@part = "+ve"
    }
  } else if(decision@stage == "4bs"){
    if(nPatOutcome0 >= p$nFinal[1]) {
      decision@analyse = TRUE
      decision@part = "-ve"
    }
  } else if(nPat > 400){
    decision@analyse = TRUE
    decision@part = "both"
  } else {
    stop("no stage found")
  }

  if(decision@stage == "1aa") {

    if(nPat0 >= p$nInt[1] & nPat1 >= p$nInt[2]) {
      decision@recruit = "none"
    } else if(nPat0 >= p$nInt[1]) {
      decision@recruit = "+ve"
    } else if(nPat1 >= p$nInt[2]){
      decision@recruit = "-ve"
    }
  } else if(decision@stage == "2ab"){

    if(nPat0 >= p$nInt[1] & nPat1 >= p$nFinal[2]) {
      decision@recruit = "none"
    } else if(nPat0 >= p$nInt[1]) {
      decision@recruit = "+ve"
    } else if(nPat1 >= p$nFinal[2]) {
      decision@recruit = "-ve"
    }
  } else if(decision@stage == "2ba"){

    if(nPat0 >= p$nFinal[1] & nPat1 >= p$nInt[2]) {
      decision@recruit = "none"
    } else if(nPat0 >= p$nFinal[1]) {
      decision@recruit = "+ve"
    } else if(nPat1 >= p$nInt[2]) {
      decision@recruit = "-ve"
    }
  } else if(decision@stage == "3bb"){
    if(nPat0 >= p$nFinal[1] & nPat1 >= p$nFinal[2]) {
      decision@recruit = "none"
    } else if(nPat0 >= p$nFinal[1]) {
      decision@recruit = "+ve"
    } else if(nPat1 >= p$nFinal[2]) {
      decision@recruit = "-ve"
    }
  } else if(decision@stage == "3fa"){
    if(nPat1 >= p$nInt[2]) {
      decision@recruit = "none"
    }
  } else if(decision@stage == "3as"){
    if(nPat0 >= p$nInt[1]) {
      decision@recruit = "none"
    }
  } else if(decision@stage == "4fb"){
    if(nPat1 >= p$nFinal[2]) {
      decision@recruit = "none"
    }
  } else if(decision@stage == "4bs"){
    if(nPat0 >= p$nFinal[1]) {
      decision@recruit = "none"
    }
  } else if(nPat > 400){
    decision@analyse = TRUE
    decision@part = "both"
  } else {
    stop("no stage found")
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

  nPat = sum(!is.na(data$outcome))
  nPat0 = sum(!is.na(data$outcome) & data$biomarker == 0)
  nPat1 = sum(!is.na(data$outcome) & data$biomarker == 1)
  # print(p$decision@stage)
  # print(p$decision@part)
  if(p$decision@stage == "1aa"){
      if(p$decision@part == "both"){

      # first interim
      # test b- and b+
      if(postProb(model$postAll,p$biomarkerH1.0, 0) < p$postProbIAll | postProb(model$postAll,p$biomarkerH0.0, 1) > p$postProb){
        if((nPat0 < p$nFinal[1] & postProb(model$post0,p$biomarkerH1.0, 0) < p$postProbIF) | (nPat0 >= p$nFinal[1] & postProb(model$post0,p$biomarkerH0.0, 1) > p$postProb)){
          stage = "3bb"
          trail = c(trail, stage)
          recruit = "both"
        } else {
          if((nPat1 < p$nFinal[2] & postProb(model$post1,p$biomarkerH1.1, 0) < p$postProbIF) | (nPat1 >= p$nFinal[2] & postProb(model$post1,p$biomarkerH0.1, 1) > p$postProb)){
            if((nPat1 >= p$nFinal[2] & postProb(model$post1,p$biomarkerH0.1, 1) > p$postProb)){
              continue = FALSE
              stopReason = "1aa success 1"
            }
            stage = "4fb"
            trail = c(trail, stage)
            recruit = "+ve"
          } else {
            continue = FALSE
            stopReason = "1aa futility"
          }
        }
      } else if((nPat1 < p$nFinal[2] & postProb(model$post1,p$biomarkerH1.1, 0) < p$postProbIF) | (nPat1 >= p$nFinal[2] & postProb(model$post1,p$biomarkerH0.1, 1) > p$postProb)){
        if((nPat1 >= p$nFinal[2] & postProb(model$post1,p$biomarkerH0.1, 1) > p$postProb)){
          continue = FALSE
          stopReason = "1aa success 1"
        }
        stage = "4fb"
        trail = c(trail, stage)
        recruit = "+ve"
      } else {
        continue = FALSE
        stopReason = "1aa futility"
      }

      } else if(p$decision@part == "-ve"){
        if((postProb(model$postAll,p$biomarkerH1.0, 0) < p$postProbIAll | postProb(model$postAll,p$biomarkerH0.0, 1) > p$postProb) & (nPat0 < p$nFinal[1] & postProb(model$post0,p$biomarkerH1.0, 0) < p$postProbIF) | (nPat0 >= p$nFinal[1] & postProb(model$post0,p$biomarkerH0.0, 1) > p$postProb)){
          stage = "2ba"
          trail = c(trail, stage)
          recruit = "both"
        } else {
          stage = "3fa"
          trail = c(trail, stage)
          recruit = "+ve"
        }
      } else if(p$decision@part == "+ve"){
        if((postProb(model$postAll,p$biomarkerH1.0, 0) < p$postProbIAll | postProb(model$postAll,p$biomarkerH0.0, 1) > p$postProb) | (nPat1 < p$nFinal[2] & postProb(model$post1,p$biomarkerH1.1, 0) < p$postProbIF) | (nPat1 >= p$nFinal[2] & postProb(model$post1,p$biomarkerH0.1, 1) > p$postProb)){
          stage = "2ab"
          trail = c(trail, stage)
          recruit = "both"
        } else {
          continue = FALSE
          stopReason = "1aa futility"
        }
      } else {
        stop("decision part not found")
      }
    } else if(p$decision@stage == "2ab"){
      if(p$decision@part == "both"){

        if((postProb(model$postAll,p$biomarkerH1.0, 0) < p$postProbIAll | postProb(model$postAll,p$biomarkerH0.0, 1) > p$postProb) & postProb(model$post0,p$biomarkerH1.0, 0) < p$postProbIF){

          if(postProb(model$post1,p$biomarkerH0.1, 1) > p$postProb){
            stage = "4bs"
            trail = c(trail, stage)
            recruit = "-ve"
          } else {
            continue = FALSE
            stopReason = "2ab futility"
          }
        } else {
          if(postProb(model$post1,p$biomarkerH0.1, 1) > p$postProb){
            continue = FALSE
            stopReason = "2ab success 1"
          } else {
            continue = FALSE
            stopReason = "2ab futility"
          }
        }
      } else if(p$decision@part == "-ve"){
        if((postProb(model$postAll,p$biomarkerH1.0, 0) < p$postProbIAll | postProb(model$postAll,p$biomarkerH0.0, 1) > p$postProb) & postProb(model$post0,p$biomarkerH1.0, 0) < p$postProbIF){
          stage = "3bb"
          trail = c(trail, stage)
          recruit = "both"
        } else {
          stage = "4fb"
          trail = c(trail, stage)
          recruit = "+ve"
        }
      } else if(p$decision@part == "+ve"){
        if(postProb(model$post1,p$biomarkerH0.1, 1) > p$postProb){
          stage = "3as"
          trail = c(trail, stage)
          recruit = "-ve"
        } else {
          continue = FALSE
          stopReason = "2ab futility"
        }
      } else {
        stop("decision part not found")
      }
    } else if(p$decision@stage == "2ba"){
      if(p$decision@part == "both"){
        if(postProb(model$postAll,p$biomarkerH0.0, 1) > p$postProb & postProb(model$post0,p$biomarkerH0.0, 1) > p$postProb){
          continue = FALSE
          stopReason = "2ba success 01"
        } else if((postProb(model$postAll,p$biomarkerH1.0, 0) < p$postProbIAll | postProb(model$postAll,p$biomarkerH0.0, 1) > p$postProb) | postProb(model$post1,p$biomarkerH1.1, 0) < p$postProbIF){
          stage = "4fb"
          trail = c(trail, stage)
          recruit = "+ve"
        } else {
          continue = FALSE
          stopReason = "2ba futility"
        }
      } else if(p$decision@part == "-ve"){
        if(postProb(model$postAll,p$biomarkerH0.0, 1) > p$postProb & postProb(model$post0,p$biomarkerH0.0, 1) > p$postProb){
          continue = FALSE
          stopReason = "4ba success 01"
        } else {
          stage = "3fa"
          trail = c(trail, stage)
          recruit = "+ve"
        }
      } else if(p$decision@part == "+ve"){
        if((postProb(model$postAll,p$biomarkerH1.0, 0) < p$postProbIAll | postProb(model$postAll,p$biomarkerH0.0, 1) > p$postProb) | postProb(model$post1,p$biomarkerH1.1, 0) < p$postProbIF){
          stage = "3bb"
          trail = c(trail, stage)
          recruit = "both"
        } else {
          continue = FALSE
          stopReason = "2ba futility"
        }
      } else {
        stop("decision part not found")
      }
    } else if(p$decision@stage == "3bb"){
      if(p$decision@part == "both"){
        if(postProb(model$postAll,p$biomarkerH0.0, 1) > p$postProb){
          if(postProb(model$post0,p$biomarkerH0.0, 1) > p$postProb){
            continue = FALSE
            stopReason = "3bb success 01"
          } else {
            if(postProb(model$post1,p$biomarkerH0.1, 1) > p$postProb){
              continue = FALSE
              stopReason = "3bb success 1"
            } else {
              continue = FALSE
              stopReason = "3bb futility"
            }
          }
        } else {
          if(postProb(model$post1,p$biomarkerH0.1, 1) > p$postProb){
            continue = FALSE
            stopReason = "3bb success 1"
          } else {
            continue = FALSE
            stopReason = "3bb futility"
          }
        }
      } else if(p$decision@part == "-ve"){
        if(postProb(model$postAll,p$biomarkerH0.0, 1) > p$postProb & postProb(model$post0,p$biomarkerH0.0, 1) > p$postProb){
          continue = FALSE
          stopReason = "3bb success 01"
        } else {
          stage = "4fb"
          trail = c(trail, stage)
          recruit = "+ve"
        }
      } else if(p$decision@part == "+ve"){
        if(postProb(model$post1,p$biomarkerH0.1, 1) > p$postProb){
          stage = "4bs"
          trail = c(trail, stage)
          recruit = "-ve"
        } else {
          continue = FALSE
          stopReason = "3bb futility"
        }
      } else {
        stop("decision part not found")
      }
    } else if(p$decision@stage == "3fa"){
      if((postProb(model$postAll,p$biomarkerH1.0, 0) < p$postProbIAll | postProb(model$postAll,p$biomarkerH0.0, 1) > p$postProb) | (nPat1 < p$nFinal[2] & postProb(model$post1,p$biomarkerH1.1, 0) < p$postProbIF) | (nPat1 >= p$nFinal[2] & postProb(model$post1,p$biomarkerH0.1, 1) > p$postProb)){
        stage = "4fb"
        trail = c(trail, stage)
        recruit = "+ve"
      } else {
        continue = FALSE
        stopReason = "3fa futility"
      }
    } else if(p$decision@stage == "3as"){
      if((postProb(model$postAll,p$biomarkerH1.0, 0) < p$postProbIAll | postProb(model$postAll,p$biomarkerH0.0, 1) > p$postProb) & (nPat0 < p$nFinal[1] & postProb(model$post0,p$biomarkerH1.0, 0) < p$postProbIF) | (nPat0 >= p$nFinal[1] & postProb(model$post0,p$biomarkerH0.0, 1) > p$postProb)){
        stage = "4bs"
        trail = c(trail, stage)
        recruit = "-ve"
      } else {
        continue = FALSE
        stopReason = "3as success 1"
      }
    } else if(p$decision@stage == "4fb"){
      if(postProb(model$post1,p$biomarkerH0.1, 1) > p$postProb){
        continue = FALSE
        stopReason = "4fb success 1"
      } else {
        continue = FALSE
        stopReason = "4fb futility"
      }
    } else if(p$decision@stage == "4bs"){
      if(postProb(model$post0,p$biomarkerH0.0, 1) > p$postProb){
        continue = FALSE
        stopReason = "4bs success 01"
      } else {
        continue = FALSE
        stopReason = "4bs success 1"
      }
  } else {
    stop("what stage is this?")
  }

  if(nPat > 400){
    continue = FALSE
  }

 # cat("\n",continue, stopReason, stage)

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
    decision = stageDecisionTwo(stage = "1aa", part = "", recruit = "both", analyse = FALSE, trail = "1aa"), # this says we want to start a stage 1a. You shouldn't need to change this.
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
  res = simTrial(setup, nSim = 1000L, parallel = TRUE)
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







