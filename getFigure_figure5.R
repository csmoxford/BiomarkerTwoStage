

inDir = paste0(outPath,"Figure 5/")
outDir = programDir

source(paste0(programDir, 'Programs/utilFunctions/utils.R'))
source(paste0(programDir,'Programs/utilFunctions/roundWZero.R'))
source(paste0(programDir, 'Programs/utilFunctions/plot_code.R'))

dfnames = c("Design","P(+)","P(b-)","P(b+)","interim stop","interim unselected","interim positive","Futility","Success (b- & b+)","Success (b+)","mean SS","99% SS","mean tested")
df = data.frame(matrix(0,nrow=0,ncol=length(dfnames)), stringsAsFactors = FALSE)
names(df) = dfnames

# Get the summary data from the saved simulations
files = list.files(inDir)
for(i in 1:length(files)){
  load(paste0(inDir,files[i]))

  stopReasons = sapply(res@sims,function(x) x$decision@stopReason)
  nrecruited = sapply(res@sims, function(x) return(sum(!is.na(x$data$outcome))))
  nseen = sapply(res@sims, function(x) return(dim(x$data)[1]))

  recommend = rep("futile",length(stopReasons))
  recommend[grepl("success 01",stopReasons)] = "group 01"
  recommend[grepl("success 1",stopReasons)] = "group 1"

  interim = rep("stop",length(stopReasons))
  interim[grepl("1b",stopReasons)] = "continue unselected"
  interim[grepl("2",stopReasons)] = "continue positive"

  res@p$design = strsplit(files[i],",")[[1]][1]
  df[i,] = c(res@p$design,res@p$biomarkerProb[2],res@p$biomarkerResp,roundWZero(mean("stop"==interim),3),roundWZero(mean("continue unselected"==interim),3),roundWZero(mean("continue positive"==interim),3),roundWZero(mean("futile"==recommend),3),roundWZero(mean("group 01"==recommend),3),roundWZero(mean("group 1"==recommend),3), roundWZero(mean(nrecruited),1), roundWZero(quantile(nrecruited,0.99)), roundWZero(mean(nseen),1))
}


#######################################################################
## Enrichment design
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


prevalences = c(0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8)

for(pr.positive in prevalences){
  method = 'opt'
  design = seq(p0,p1,alpha1,power1,alpha0,power0,n.max,u,method)
  df = rbind(df, enrichGetPropertiesTable(design, "Enrichment first (Optimal)",pr.positive))
  method = 'min'
  design = seq(p0,p1,alpha1,power1,alpha0,power0,n.max,u,method)
  df = rbind(df, enrichGetPropertiesTable(design, "Enrichment first (Minmax)",pr.positive))
}


#######################################################################
## Parrallel design

##############################################################################
# Independent parallel
source(paste0(programDir, 'Programs/Independent_parallel.R'))

p0=0.1
p1=0.3
n1 = 17
n = 34

eta=c(0.98,0.95)
zeta=c(NA,0.90)
prior.a = 0.2/4
prior.b = 0.8/4

for(pr.positive in prevalences){
  df = rbind(df,getIndParallelProperties("IndParallel",pr.positive, p0,p1, res, n1, n, eta, zeta, prior.a,prior.b))
}


st = df

save(df, file = paste0(outPath,"Figure5_summary.rData"))
##################################################
load(paste0(outPath,"Figure5_summary.rData"))
the_design_names = unique(df$Design)

the_design_names = c( "IndParallel","two_betaBinom_tandem","two_betaBinom_balanced", "Enrichment first (Optimal)", "Enrichment first (Minmax)")

print_names = c("Independent Parallel","Tandem two stage","Parallel two-stage balanced", "Sequential Enrichment (optimal)", "Sequential Enrichment (minmax)")

# colours
cols = cols = brewer.pal(6, "Accent")[c(1,2,3,5,6)]
df$cols = sapply(df$Design, function(x) cols[which( x == the_design_names)][1])

# create subset for each hypothesis scenario
df1 = df[df$`P(b-)`== 0.1 & df$`P(b+)` == 0.1,]
df2 = df[df$`P(b-)`== 0.1 & df$`P(b+)` == 0.3,]
df3 = df[df$`P(b-)`== 0.3 & df$`P(b+)` == 0.3,]

# pdf(paste0(outDir,"Figure 5_colorCode.pdf"), width = 6, height = 6)
par(mar=c(0,0,0,0))
plot(rep(0.5,length(the_design_names)),1:length(the_design_names), xlim=c(0,5), col = cols, pch = 16, cex = 3, axes = FALSE)
text(0.8,1:length(the_design_names), labels = print_names, pos = 4)
box()
# dev.off()


# save to eps rather than pdf
# pdf(paste0(outDir,"Figure 5.pdf"), width = 14, height = 14/3)
setEPS(width = 14, height = 14/3, pointsize = 12)
postscript(paste0(outDir,"Figure 5.eps"))

split.screen(rbind(c(0.03,0.80,0.88,1),c(0.03,0.80,0.05,0.88),c(0.80,1,0.05,0.88)))
screen(1)
sc = split.screen(c(1,3))
screen(sc[1])
par(mar=c(0,3,0,0.1))
blankPlot()
text(0.5,0.5, label = expression(P(R^{"-"}) ~ "=" ~ P(R^{"+"}) ~ "= 0.1"), cex = 1.5, adj = 0.5)
screen(sc[2])
par(mar=c(0,3,0,0.1))
blankPlot()
text(0.5,0.5, label = expression(P(R^{"-"}) ~ "= 0.1, " ~ P(R^{"+"}) ~ "= 0.3"), cex = 1.5, adj = 0.5)
screen(sc[3])
par(mar=c(0,3,0,0.1))
blankPlot()
text(0.5,0.5, label = expression(P(R^{"-"}) ~ "=" ~ P(R^{"+"}) ~ "= 0.3"), cex = 1.5, adj = 0.5)

sc = split.screen(c(1,3),2)

ylim = range(pretty(as.numeric(df$`mean tested`)))

screen(sc[1])
plotN(df1, ylim = ylim, plotMean = FALSE, ylab = TRUE, xlab = TRUE)
screen(sc[2])
plotN(df2, ylim = ylim, plotMean = FALSE, ylab = FALSE, xlab = TRUE)
screen(sc[3])
plotN(df3, ylim = ylim, plotMean = FALSE, ylab = FALSE, xlab = TRUE)


screen(3)
par(mar=c(0,0,0,0))
plot(rep(0.3,length(the_design_names)),length(the_design_names):1, xlim=c(0,5), col = cols, pch = 16, cex = 2.5, axes = FALSE, ylim = c(-1,6))
text(0.6,length(the_design_names):1, labels = print_names, pos = 4, cex = 0.95)


close.screen(all.screens = TRUE)

dev.off()
