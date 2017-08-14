
## Run the simulations for table 1 of the paper.

outDir = paste0(outPath,"Table 1/")


n = 10000L

###############################################################################
design = "two_betaBinom_tandem"
source(paste0(programDir, 'Programs/', design,'.R'))

setup@seed = 5356L

setup@p$biomarkerResp = c(0.1,0.1)
setup@p$postProb = 0.95

res = simTrial(setup,nSim = n, parallel = TRUE)
save(res,file = paste0(outDir,design,",",paste0(setup@p$biomarkerResp,collapse = ","),".Rdata"))

setup@p$biomarkerResp = c(0.1,0.3)

res = simTrial(setup,nSim = n, parallel = TRUE)
save(res,file = paste0(outDir,design,",",paste0(setup@p$biomarkerResp,collapse = ","),".Rdata"))

setup@p$biomarkerResp = c(0.3,0.3)

res = simTrial(setup,nSim = n, parallel = TRUE)
save(res,file = paste0(outDir,design,",",paste0(setup@p$biomarkerResp,collapse = ","),".Rdata"))


###############################################################################
design = "two_betaBinom_bestFirst"
source(paste0(programDir, 'Programs/', design,'.R'))

setup@seed = 5356L

setup@p$biomarkerResp = c(0.1,0.1)
setup@p$postProb = 0.95

res = simTrial(setup,nSim = n, parallel = TRUE)
save(res,file = paste0(outDir,design,",",paste0(setup@p$biomarkerResp,collapse = ","),".Rdata"))

setup@p$biomarkerResp = c(0.1,0.3)

res = simTrial(setup,nSim = n, parallel = TRUE)
save(res,file = paste0(outDir,design,",",paste0(setup@p$biomarkerResp,collapse = ","),".Rdata"))

setup@p$biomarkerResp = c(0.3,0.3)

res = simTrial(setup,nSim = n, parallel = TRUE)
save(res,file = paste0(outDir,design,",",paste0(setup@p$biomarkerResp,collapse = ","),".Rdata"))

###############################################################################
design = "two_betaBinom_worstFirst"
source(paste0(programDir, 'Programs/', design,'.R'))

setup@seed = 5356L

setup@p$biomarkerResp = c(0.1,0.1)
setup@p$postProb = 0.95

res = simTrial(setup,nSim = n, parallel = TRUE)
save(res,file = paste0(outDir,design,",",paste0(setup@p$biomarkerResp,collapse = ","),".Rdata"))

setup@p$biomarkerResp = c(0.1,0.3)

res = simTrial(setup,nSim = n, parallel = TRUE)
save(res,file = paste0(outDir,design,",",paste0(setup@p$biomarkerResp,collapse = ","),".Rdata"))

setup@p$biomarkerResp = c(0.3,0.3)

res = simTrial(setup,nSim = n, parallel = TRUE)
save(res,file = paste0(outDir,design,",",paste0(setup@p$biomarkerResp,collapse = ","),".Rdata"))


###############################################################################
design = "two_betaBinom_balanced"
source(paste0(programDir, 'Programs/', design,'.R'))

source(paste0(programDir, 'Programs/staggered interim/two_betaBinom_balanced_flexible_pauseAll.R'))

setup@seed = 5356L
setup@p$stopRule = "pause"
setup@p$biomarkerResp = c(0.1,0.1)
setup@p$postProb = 0.95

res = simTrial(setup,nSim = n, parallel = TRUE)
save(res,file = paste0(outDir,design,",",paste0(setup@p$biomarkerResp,collapse = ","),".Rdata"))

setup@p$biomarkerResp = c(0.1,0.3)

res = simTrial(setup,nSim = n, parallel = TRUE)
save(res,file = paste0(outDir,design,",",paste0(setup@p$biomarkerResp,collapse = ","),".Rdata"))

setup@p$biomarkerResp = c(0.3,0.3)

res = simTrial(setup,nSim = n, parallel = TRUE)
save(res,file = paste0(outDir,design,",",paste0(setup@p$biomarkerResp,collapse = ","),".Rdata"))

###############################################################################
design = "two_logisticRestricted_bestFirst"
source(paste0(programDir, 'Programs/', design,'.R'))

setup@seed = 5356L

setup@p$biomarkerResp = c(0.1,0.1)
setup@p$postProb = 0.95

res = simTrial(setup,nSim = n, parallel = TRUE)
save(res,file = paste0(outDir,design,",",paste0(setup@p$biomarkerResp,collapse = ","),".Rdata"))

setup@p$biomarkerResp = c(0.1,0.3)

res = simTrial(setup,nSim = n, parallel = TRUE)
save(res,file = paste0(outDir,design,",",paste0(setup@p$biomarkerResp,collapse = ","),".Rdata"))

setup@p$biomarkerResp = c(0.3,0.3)

res = simTrial(setup,nSim = n, parallel = TRUE)
save(res,file = paste0(outDir,design,",",paste0(setup@p$biomarkerResp,collapse = ","),".Rdata"))
