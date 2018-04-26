
# Run simulations for figure 2

outDir = paste0(outPath,"Figure 5/")

n = 10000L


prevalences = c(0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8)
###############################################################################
design = "two_betaBinom_tandem"
source(paste0(programDir, 'Programs/', design,'.R'))

setup@seed = 5356L

setup@p$postProb = 0.95

for(pr.positive in prevalences){

  setup@p$biomarkerProb = c(1-pr.positive,pr.positive)

  setup@p$biomarkerResp = c(0.1,0.1)

  res = simTrial(setup,nSim = n, parallel = TRUE)
  save(res,file = paste0(outDir,design,",",paste0(setup@p$biomarkerResp,collapse = ","),",",pr.positive,".Rdata"))

  setup@p$biomarkerResp = c(0.1,0.3)

  res = simTrial(setup,nSim = n, parallel = TRUE)
  save(res,file = paste0(outDir,design,",",paste0(setup@p$biomarkerResp,collapse = ","),",",pr.positive,".Rdata"))

  setup@p$biomarkerResp = c(0.3,0.3)

  res = simTrial(setup,nSim = n, parallel = TRUE)
  save(res,file = paste0(outDir,design,",",paste0(setup@p$biomarkerResp,collapse = ","),",",pr.positive,".Rdata"))

}


###############################################################################
design = "two_betaBinom_bestFirst"
source(paste0(programDir, 'Programs/', design,'.R'))

setup@seed = 5356L

setup@p$postProb = 0.95

for(pr.positive in prevalences){

  setup@p$biomarkerProb = c(1-pr.positive,pr.positive)

  setup@p$biomarkerResp = c(0.1,0.1)

  res = simTrial(setup,nSim = n, parallel = TRUE)
  save(res,file = paste0(outDir,design,",",paste0(setup@p$biomarkerResp,collapse = ","),",",pr.positive ,".Rdata"))

  setup@p$biomarkerResp = c(0.1,0.3)

  res = simTrial(setup,nSim = n, parallel = TRUE)
  save(res,file = paste0(outDir,design,",",paste0(setup@p$biomarkerResp,collapse = ","),",",pr.positive,".Rdata"))

  setup@p$biomarkerResp = c(0.3,0.3)

  res = simTrial(setup,nSim = n, parallel = TRUE)
  save(res,file = paste0(outDir,design,",",paste0(setup@p$biomarkerResp,collapse = ","),",",pr.positive,".Rdata"))

}

###############################################################################
design = "two_betaBinom_worstFirst"
source(paste0(programDir, 'Programs/', design,'.R'))

setup@seed = 5356L

setup@p$postProb = 0.95

for(pr.positive in prevalences){

  setup@p$biomarkerProb = c(1-pr.positive,pr.positive)

  setup@p$biomarkerResp = c(0.1,0.1)

  res = simTrial(setup,nSim = n, parallel = TRUE)
  save(res,file = paste0(outDir,design,",",paste0(setup@p$biomarkerResp,collapse = ","),",",pr.positive,".Rdata"))

  setup@p$biomarkerResp = c(0.1,0.3)

  res = simTrial(setup,nSim = n, parallel = TRUE)
  save(res,file = paste0(outDir,design,",",paste0(setup@p$biomarkerResp,collapse = ","),",",pr.positive,".Rdata"))

  setup@p$biomarkerResp = c(0.3,0.3)

  res = simTrial(setup,nSim = n, parallel = TRUE)
  save(res,file = paste0(outDir,design,",",paste0(setup@p$biomarkerResp,collapse = ","),",",pr.positive,".Rdata"))

}

###############################################################################
design = "two_betaBinom_balanced"
source(paste0(programDir, 'Programs/', design,'.R'))

setup@seed = 5356L

setup@p$postProb = 0.95
for(pr.positive in prevalences){

  setup@p$biomarkerProb = c(1-pr.positive,pr.positive)

  setup@p$biomarkerResp = c(0.1,0.1)

  res = simTrial(setup,nSim = n, parallel = TRUE)
  save(res,file = paste0(outDir,design,",",paste0(setup@p$biomarkerResp,collapse = ","),",",pr.positive,".Rdata"))

  setup@p$biomarkerResp = c(0.1,0.3)

  res = simTrial(setup,nSim = n, parallel = TRUE)
  save(res,file = paste0(outDir,design,",",paste0(setup@p$biomarkerResp,collapse = ","),",",pr.positive,".Rdata"))

  setup@p$biomarkerResp = c(0.3,0.3)

  res = simTrial(setup,nSim = n, parallel = TRUE)
  save(res,file = paste0(outDir,design,",",paste0(setup@p$biomarkerResp,collapse = ","),",",pr.positive,".Rdata"))

  }

###############################################################################
design = "two_logisticRestricted_bestFirst"
source(paste0(programDir, 'Programs/', design,'.R'))

setup@seed = 5356L

setup@p$postProb = 0.95

for(pr.positive in prevalences){

  setup@p$biomarkerProb = c(1-pr.positive,pr.positive)

  setup@p$biomarkerResp = c(0.1,0.1)

  res = simTrial(setup,nSim = n, parallel = TRUE)
  save(res,file = paste0(outDir,design,",",paste0(setup@p$biomarkerResp,collapse = ","),",",pr.positive,".Rdata"))

  setup@p$biomarkerResp = c(0.1,0.3)

  res = simTrial(setup,nSim = n, parallel = TRUE)
  save(res,file = paste0(outDir,design,",",paste0(setup@p$biomarkerResp,collapse = ","),",",pr.positive,".Rdata"))

  setup@p$biomarkerResp = c(0.3,0.3)

  res = simTrial(setup,nSim = n, parallel = TRUE)
  save(res,file = paste0(outDir,design,",",paste0(setup@p$biomarkerResp,collapse = ","),",",pr.positive,".Rdata"))

}
