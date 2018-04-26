###############################################################
# Stratified medicine: Non-randomised two-stage designs for molecularly targeted treatments
# P. Dutton, J. Holmes
###############################################################
# The trialSim package provides a framework for running both simple and more complext trial simulations with individual patient data to assist in optimising and obtaining the properties of different trial designs with the eventual goal to choose one for a clinical trial or experiment. The package is designed to make everything reproducible and allows easy parallelisation of individual simulations.
# The classes are generated in each of .R files in the programs folder with the exception of the independent parallel and the Enrichment design which does not use the trialSim package to generate the running properties of the designs.
# The code on this page should run provided the programDir points to the folder containing this file and the outPath points to a valid folder. Note that all simulations from trialSim are saved and as a result the saved simulation are around 7GB with a long runtime to generate figures 5 and 6. Not that figure 5 is figure 6 in the paper and figure 6 is figure 8 in the paper
###############################################################
library(trialSim) # devtools::install_github("finite2/trialSim")
library(prettyTables) # devtools::install_github("csmoxford/prettyTables")
library(clinfun)
library(knitr)
library(rmarkdown)
library(R2jags)
library(stringr)
library(EurosarcBayes)
library(RColorBrewer)

# Top level of this directory
programDir = "I:/Data/Peter_Dutton/Papers/wee1Paper/BiomarkerTwoStage/"
# save path if difference from programDir
outPath = "H:/wee1 paper/Simulations_15May2017/"
###############################################################
## Run exact calculations for table 1 of the paper

# The code underlying table 1 is in:
# getData_table1_18Dec2017.R
## This calls a number of routines:
# TandemTwoStage_GetExactProperties.R
# ParallelTwoStage_GetExactProperties.R
# enrich_GetExactPropertiesTable.R

# Compile table 1 to a word document
render(
  input = paste0(programDir,"getTable_table1_02Feb2018.Rmd"),
  output_file = paste0(programDir,"Table 1.docx")
)

###############################################################
## Run the simulations for figure 6 of the paper
source(paste0(programDir,"getFigure_figure6_20Dec2017"), echo = TRUE)

###############################################################
## Run the simulations for figure 8 of the paper
source(paste0(programDir,"run_figure8.R"), echo = TRUE)

# Compile summary data and produce figure 8
source(paste0(programDir,"getFigure_figure8.R"), echo = TRUE)

###############################################################

