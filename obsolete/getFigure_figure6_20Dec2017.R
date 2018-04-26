   

load(paste0(programDir, "BalancedStaggeredData_02Feb2018.rData"))
df = rbind(df, balancedData)
the_design_names = unique(df$Design)

the_design_names = c( "IndParallel","tandemTwoStage", "worstFirst","Balanced Staggered", "Sequential Enrichment")

print_names = c("Independent Parallel","Tandem two stage", "Parallel two-stage negative-first","Parallel two-stage balanced", "Sequential Enrichment")

# colours
cols = cols = brewer.pal(7, "Accent")[c(1,2,3,5,6)]
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



source(paste0(programDir,"Programs/utilFunctions/plot_code.R"))
# save to eps rather than pdf
# pdf(paste0(outDir,"Figure 5.pdf"), width = 14, height = 14/3)
setEPS(width = 14, height = 14/3, pointsize = 12)
postscript(paste0(programDir,"Figure 5.eps"))

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
