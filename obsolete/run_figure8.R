
n = 10000L

outDir = paste0(outPath,"/Figure 6/")

source(paste0(programDir,"Programs/staggered interim/simulationSetup.R"))

###############################################################################

pauseAll = FALSE
for(probPos in c(0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8)){

  source(paste0(programDir, 'Programs/staggered interim/two_betaBinom_balanced_flexible_pauseAll.R'))
  ## pauseOption: No staggering of interim analysis wait to analyse
  setup@seed = 5356L

  setup@p$biomarkerProb = c(1-probPos, probPos)

  design = paste0("betabinom_balanced_all_nostagger",probPos)

  setup@p$nWait = 35

  setup@p$biomarkerResp = c(0.1,0.1)

  res = simTrial(setup,nSim = n, parallel = TRUE)
  save(res,file = paste0(outDir, design,",",paste0(setup@p$biomarkerResp,collapse = ","),".Rdata"))

  setup@p$biomarkerResp = c(0.1,0.3)

  res = simTrial(setup,nSim = n, parallel = TRUE)
  save(res,file = paste0(outDir, design,",",paste0(setup@p$biomarkerResp,collapse = ","),".Rdata"))

  setup@p$biomarkerResp = c(0.3,0.3)

  res = simTrial(setup,nSim = n, parallel = TRUE)
  save(res,file = paste0(outDir, design,",",paste0(setup@p$biomarkerResp,collapse = ","),".Rdata"))

}


for(probPos in c(0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8)){

  source(paste0(programDir, 'Programs/staggered interim/two_betaBinom_balanced_flexible_pauseAnalysis.R'))

  setup@seed = 5356L
  setup@p$nWait = 3

  ## pauseOption: pause analysis group only
  setup@p$biomarkerProb = c(1-probPos, probPos)

  design = paste0("betabinom_balanced_all_analyse_3wait",probPos)

  setup@p$biomarkerResp = c(0.1,0.1)

  res = simTrial(setup,nSim = n, parallel = TRUE)
  save(res,file = paste0(outDir, design,",",paste0(setup@p$biomarkerResp,collapse = ","),".Rdata"))

  setup@p$biomarkerResp = c(0.1,0.3)

  res = simTrial(setup,nSim = n, parallel = TRUE)
  save(res,file = paste0(outDir, design,",",paste0(setup@p$biomarkerResp,collapse = ","),".Rdata"))

  setup@p$biomarkerResp = c(0.3,0.3)

  res = simTrial(setup,nSim = n, parallel = TRUE)
  save(res,file = paste0(outDir, design,",",paste0(setup@p$biomarkerResp,collapse = ","),".Rdata"))

}


for(probPos in c(0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8)){

  source(paste0(programDir, 'Programs/staggered interim/two_betaBinom_balanced_flexible_pauseAnalysis.R'))

  setup@seed = 5356L
  setup@p$nWait = 0

  ## pauseOption: pause analysis group only
  setup@p$biomarkerProb = c(1-probPos, probPos)

  design = paste0("betabinom_balanced_all_analyse_nowait",probPos)

  setup@p$biomarkerResp = c(0.1,0.1)

  res = simTrial(setup,nSim = n, parallel = TRUE)
  save(res,file = paste0(outDir, design,",",paste0(setup@p$biomarkerResp,collapse = ","),".Rdata"))

  setup@p$biomarkerResp = c(0.1,0.3)

  res = simTrial(setup,nSim = n, parallel = TRUE)
  save(res,file = paste0(outDir, design,",",paste0(setup@p$biomarkerResp,collapse = ","),".Rdata"))

  setup@p$biomarkerResp = c(0.3,0.3)

  res = simTrial(setup,nSim = n, parallel = TRUE)
  save(res,file = paste0(outDir, design,",",paste0(setup@p$biomarkerResp,collapse = ","),".Rdata"))

}

for(probPos in c(0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8)){
  source(paste0(programDir, 'Programs/staggered interim/two_betaBinom_balanced_flexible_FinalOnly.R'))
  ## pauseOption: pause for final anlysis only
  setup@p$biomarkerProb = c(1-probPos, probPos)

  setup@seed = 5356L

  design = paste0("betabinom_balanced_all_final",probPos)

  setup@p$biomarkerResp = c(0.1,0.1)

  res = simTrial(setup,nSim = n, parallel = TRUE)
  save(res,file = paste0(outDir, design,",",paste0(setup@p$biomarkerResp,collapse = ","),".Rdata"))

  setup@p$biomarkerResp = c(0.1,0.3)

  res = simTrial(setup,nSim = n, parallel = TRUE)
  save(res,file = paste0(outDir, design,",",paste0(setup@p$biomarkerResp,collapse = ","),".Rdata"))

  setup@p$biomarkerResp = c(0.3,0.3)

  res = simTrial(setup,nSim = n, parallel = TRUE)
  save(res,file = paste0(outDir, design,",",paste0(setup@p$biomarkerResp,collapse = ","),".Rdata"))

}
