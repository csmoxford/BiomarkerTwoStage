
source(paste0(programDir, 'Programs/Enrichment_Design_Zang_Yuan_2016.R'))
source(paste0(programDir, 'Programs/utilFunctions/enrichGetPropertiesTable.R'))

p0 = 0.1  ## p0: unacceptable response rate
p1 = 0.3  ## p1: desirable response rate
alpha1 = 0.1  ## alpha1: type I error for marker-positive subgroup
power1 = 0.925  ## power1: power for marker-positive subgroup
alpha0 = 0.1  ## alpha0: type I error for marker-negative subgroup
power0 = 0.85  ## power0: power for marker-negative subgroup
n.max = 100  ## n.max: maximum sample size for each marker subgroup
u = 0.3  ## u: the upper bound of response rate for marker positive subgroup
method = 'opt'  ## method: "opt" refers to the optimal design; "min" refers to the minmax design


pr.positive = 0.5

design = seq(p0,p1,alpha1,power1,alpha0,power0,n.max,u,method)
dta = enrichGetPropertiesTable(design, "Enrichment first (Optimal)",pr.positive)
method = 'min'
design = seq(p0,p1,alpha1,power1,alpha0,power0,n.max,u,method)
dta1 = enrichGetPropertiesTable(design, "Enrichment first (Minmax)",pr.positive)














