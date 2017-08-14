
inDir = paste0(outPath,"Figure 6/")
outDir = programDir

source(paste0(programDir, 'Programs/utilFunctions/utils.R'))
source(paste0(programDir,'Programs/utilFunctions/roundWZero.R'))
source(paste0(programDir, 'Programs/utilFunctions/plot_code.R'))

## full summary
dfnames = c("Design","P(+)","P(b-)","P(b+)","Futility","Success (b- & b+)","Success (b+)","mean SS","99% SS","mean tested","mean Time")
df = data.frame(matrix(0,nrow=0,ncol=length(dfnames)), stringsAsFactors = FALSE)
names(df) = dfnames



files = list.files(inDir)
# Get the summary data from the saved simulations
for(i in 1:length(files)){

  load(paste0(inDir,files[i]))

  stopReasons = sapply(res@sims,function(x) x$decision@stopReason)
  nrecruited = sapply(res@sims, function(x) return(sum(!is.na(x$data$outcome))))
  nseen = sapply(res@sims, function(x) return(dim(x$data)[1]))

  time = sapply(res@sims, function(x) max(x$data$arrivalTime[!is.na(x$data$outcome)]))

  recommend = rep("futile",length(stopReasons))
  recommend[grepl("success 01",stopReasons)] = "group 01"
  recommend[grepl("success 1",stopReasons)] = "group 1"

  interim = rep("stop",length(stopReasons))
  interim[grepl("1b",stopReasons)] = "continue unselected"
  interim[grepl("2",stopReasons)] = "continue positive"

  res@p$design = strsplit(files[i],",")[[1]][1]
  df[i,] = c(res@p$design,res@p$biomarkerProb[2],res@p$biomarkerResp,roundWZero(mean("futile"==recommend),3),roundWZero(mean("group 01"==recommend),3),roundWZero(mean("group 1"==recommend),3), roundWZero(mean(nrecruited),1), roundWZero(quantile(nrecruited,0.99)), roundWZero(mean(nseen),1),roundWZero(mean(time),1))
}

save(df,file = paste0(outDir,"Figure6_data.rData"))
##############################################################
load(paste0(outDir,"/Figure6_data.rData"))

library(RColorBrewer)
table(df$Design)


the_design_names = c("betabinom_balanced_all_nostagger", "betabinom_balanced_all_analyse_threewait", "betabinom_balanced_all_analyse_nowait", "betabinom_balanced_all_final")
labels = c("No staggered interim", "Staggered, with pausing, wait for three patients", "Staggered, with pausing", "Staggered without pausing")

# colours
cols = brewer.pal(6, "Accent")[c(1,2,5,6)]

df$Design = gsub("3analyse","threeanalyse",df$Design)
df$Design = gsub("3wait","threewait",df$Design)
df$Design = gsub("[0-9.]","",df$Design)
df$cols = sapply(df$Design, function(x) cols[which( x == the_design_names)])

# create subset for each hypothesis scenario
df1 = df[df$`P(b-)`== 0.1 & df$`P(b+)` == 0.1,]
df2 = df[df$`P(b-)`== 0.1 & df$`P(b+)` == 0.3,]
df3 = df[df$`P(b-)`== 0.3 & df$`P(b+)` == 0.3,]


# pdf(paste0(outDir,"Figure 6_colorCode.pdf"), width = 6, height = 6)
par(mar=c(0,0,0,0))
plot(rep(0.2,length(labels)),1:length(labels), col = cols, pch = 16, cex = 4, axes = FALSE,xlim=c(0,1))
text(0.6,1:length(labels), labels = labels, cex = 1.5)
box()
# dev.off()

# create subset for each hypothesis scenario

# save to eps rather than pdf
# pdf(paste0(outDir,"Figure 6.pdf"), width = 14, height = 14)
setEPS(width = 14, height = 14, pointsize = 12)
postscript(paste0(outDir,"Figure 6.eps"))


split.screen(rbind(c(0.05,0.98,0.96,1),c(0.05,0.98,0.05,0.96)))
screen(1)
split.screen(c(1,3))
screen(3)
par(mar=c(0,2,0,0))
plot(0,0,xlim=c(0,1),ylim=c(0,1),type = "n", axes = FALSE, xlab = "", ylab = "")
text(1/2,0.5, label = expression(P(R^{"-"}) ~ "=" ~ P(R^{"+"}) ~ "= 0.1"), cex = 1.5, adj = 0.5)
screen(4)
par(mar=c(0,2,0,0))
plot(0,0,xlim=c(0,1),ylim=c(0,1),type = "n", axes = FALSE, xlab = "", ylab = "")
text(1/2,0.5, label = expression(P(R^{"-"}) ~ "= 0.1, " ~ P(R^{"+"}) ~ "= 0.3"), cex = 1.5, adj = 0.5)
screen(5)
par(mar=c(0,2,0,0))
plot(0,0,xlim=c(0,1),ylim=c(0,1),type = "n", axes = FALSE, xlab = "", ylab = "")
text(1/2,0.5, label = expression(P(R^{"-"}) ~ "=" ~ P(R^{"+"}) ~ "= 0.3"), cex = 1.5, adj = 0.5)

sc = split.screen(c(3,3),2)

ylim = c(0.82,0.99)

screen(sc[1])
plotTruePositive(df1, "Futility", ylab = TRUE, ylim = ylim)
screen(sc[2])
plotTruePositive(df2, "Success (b+)",ylab = FALSE, ylim = ylim)
screen(sc[3])
plotTruePositive(df3, "Success (b- & b+)",ylab = FALSE, ylim = ylim)

ylim = range(pretty(as.numeric(df[,"mean SS"])))

screen(sc[4])
plotN(df1,tested = FALSE, ylim = ylim)
screen(sc[5])
plotN(df2,tested = FALSE,ylab = FALSE, ylim = ylim)
screen(sc[6])
plotN(df3,tested = FALSE,ylab = FALSE, ylim = ylim)
legend("bottomright",pch = 16,pt.cex = 2, col = cols[1:4], legend = labels, cex = 0.9)

ylim = range(pretty(as.numeric(df[,"mean Time"])))

screen(sc[7])
plotTime(df1, ylim = ylim)
screen(sc[8])
plotTime(df2,ylab = FALSE, ylim = ylim)
screen(sc[9])
plotTime(df3,ylab = FALSE, ylim = ylim)

close.screen(all.screens = TRUE)

dev.off()
