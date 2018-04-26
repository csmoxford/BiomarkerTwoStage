
inDir = paste0(outPath,"Figure 6/")

source(paste0(programDir, 'Programs/utilFunctions/utils.R'))
source(paste0(programDir,'Programs/utilFunctions/roundWZero.R'))
source(paste0(programDir, 'Programs/utilFunctions/plot_code.R'))

## full summary
dfnames = c("Design","P(+)","P(b-)","P(b+)","Futility","Success (b- & b+)","Success (b+)","mean SS","99% SS","mean tested","mean Time")
df = data.frame(matrix(0,nrow=0,ncol=length(dfnames)), stringsAsFactors = FALSE)
names(df) = dfnames

## get data

load(paste0(programDir, "Figure5_summary.rData"))


load(paste0(programDir, "BalancedStaggeredData_02Feb2018.rData"))
df = rbind(df,balancedData)
##



library(RColorBrewer)
table(df$Design)


the_design_names = c("IndParallel", "worstFirst", "balanced", "Balanced Staggered", "Sequential Enrichment")
labels = c("No staggered interim", "Staggered, with pausing, wait for three patients", "Staggered, with pausing", "Staggered without pausing")

labels = c("Independent Parallel", "Parallel two-stage negative-first","Parallel two-stage balanced (wait)", "Parallel two-stage balanced", "Sequential Enrichment")

df = df[df$Design %in% the_design_names, ]

# colours
cols = brewer.pal(7, "Accent")[c(1,7,5,6)]

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
setEPS(width = 14, height = 14/3, pointsize = 12)
postscript(paste0(outDir,"Figure 6.eps"))


split.screen(rbind(c(0.03,0.80,0.88,1),c(0.03,0.80,0.05,0.88),c(0.80,1,0.05,0.88)))
screen(1)
sc = split.screen(c(1,3))
screen(sc[1])
par(mar=c(0,2,0,0))
plot(0,0,xlim=c(0,1),ylim=c(0,1),type = "n", axes = FALSE, xlab = "", ylab = "")
text(1/2,0.5, label = expression(P(R^{"-"}) ~ "=" ~ P(R^{"+"}) ~ "= 0.1"), cex = 1.5, adj = 0.5)
screen(sc[2])
par(mar=c(0,2,0,0))
plot(0,0,xlim=c(0,1),ylim=c(0,1),type = "n", axes = FALSE, xlab = "", ylab = "")
text(1/2,0.5, label = expression(P(R^{"-"}) ~ "= 0.1, " ~ P(R^{"+"}) ~ "= 0.3"), cex = 1.5, adj = 0.5)
screen(sc[3])
par(mar=c(0,2,0,0))
plot(0,0,xlim=c(0,1),ylim=c(0,1),type = "n", axes = FALSE, xlab = "", ylab = "")
text(1/2,0.5, label = expression(P(R^{"-"}) ~ "=" ~ P(R^{"+"}) ~ "= 0.3"), cex = 1.5, adj = 0.5)

sc = split.screen(c(1,3),2)

ylim = c(0.83,0.97)

screen(sc[1])
plotTruePositive(df1, "Futility", ylab = TRUE, ylim = ylim, lineIfEqual = TRUE)
screen(sc[2])
plotTruePositive(df2, "Success (b+)",ylab = FALSE, ylim = ylim, lineIfEqual = TRUE)
screen(sc[3])
plotTruePositive(df3, "Success (b- & b+)",ylab = FALSE, ylim = ylim, lineIfEqual = TRUE)

screen(3)
par(mar=c(0,0,0,0))
plot(rep(0.3,length(the_design_names)),length(the_design_names):1, xlim=c(0,5), col = cols, pch = 16, cex = 2.5, axes = FALSE, ylim = c(-1,6))
text(0.6,length(the_design_names):1, labels = print_names, pos = 4, cex = 0.95)



close.screen(all.screens = TRUE)

dev.off()
