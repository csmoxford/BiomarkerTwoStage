
inDir = paste0(outPath,"Table 1/")
outDir = outPath

source(paste0(programDir, 'Programs/utilFunctions/utils.R'))
source(paste0(programDir,'Programs/utilFunctions/roundWZero.R'))

## full summary
dfnames = c("Design","P(+)","P(b-)","P(b+)","interim stop","interim unselected","interim positive","Futility","Success (b- & b+)","Success (b+)","mean SS","99% SS","mean tested")
df = data.frame(matrix(0,nrow=0,ncol=length(dfnames)), stringsAsFactors = FALSE)
names(df) = dfnames


files = list.files(inDir)
#
for(i in 1:length(files)){
  load(paste0(inDir,files[i]))

  stopReasons = sapply(res@sims,function(x) x$decision@stopReason)
  nrecruited = sapply(res@sims, function(x) return(sum(!is.na(x$data$outcome))))
  nseen = sapply(res@sims, function(x) return(sum(!is.na(x$data$biomarker), na.rm = TRUE) - 1)) # -1 for last row.

  recommend = rep("futile",length(stopReasons))
  recommend[grepl("success 01",stopReasons)] = "group 01"
  recommend[grepl("success 1",stopReasons)] = "group 1"

  interim = rep("stop",length(stopReasons))
  interim[grepl("1b",stopReasons)] = "continue unselected"
  interim[grepl("2",stopReasons)] = "continue positive"

  res@p$design = strsplit(files[i],",")[[1]][1]
  df[i,] = c(res@p$design,res@p$biomarkerProb[2],res@p$biomarkerResp,roundWZero(mean("stop"==interim),3),roundWZero(mean("continue unselected"==interim),3),roundWZero(mean("continue positive"==interim),3),roundWZero(mean("futile"==recommend),3),roundWZero(mean("group 01"==recommend),3),roundWZero(mean("group 1"==recommend),3), roundWZero(mean(nrecruited),1), roundWZero(quantile(nrecruited,0.99)), roundWZero(mean(nseen),1))
}

##############################################################################
# Independent parallel
source(paste0(programDir,"run_table1_Independent_parallel.R"))

p0=0.1
p1=0.3
pr.positive = 0.5
n1 = 17
n = 34

eta=c(0.98,0.95)
zeta=c(NA,0.90)
prior.a = 0.2/4
prior.b = 0.8/4

df = rbind(df,getIndParallelProperties("IndParallel",pr.positive, p0,p1, res, n1, n, eta, zeta, prior.a,prior.b))

##############################################################################
# Sequential enrichment
source(paste0(programDir,"run_table1_Zang_Yuan.R"))

df = rbind(df,dta)
df = rbind(df,dta1)

df

save(df,file = paste0(outDir, "table1_data.Rdata"))






