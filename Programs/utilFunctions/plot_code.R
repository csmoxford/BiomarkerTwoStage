
blankPlot = function (xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i") {
  plot(0, 0, axes = FALSE, xlim = xlim, ylim = ylim, xlab = "",
       ylab = "", xaxs = xaxs, yaxs = yaxs, type = "n")
}

plotN = function(df1, tested = TRUE, plotMean = TRUE, ylab = TRUE, xlab = FALSE, ylim){

  par(mar=c(3.2,3,0.1,0.1))
  if(tested){
    plot(c(df1$"P(+)", df1$"P(+)"),c(df1$"mean SS", df1[,"mean tested"]),col=0, axes = FALSE, xlab = "", ylab = "", ylim = ylim)
  } else {
    plot(df1$"P(+)",df1$"mean SS",col=0, axes = FALSE, xlab = "", ylab = "", ylim = ylim)
  }
  for(design in the_design_names){
    if(plotMean){
    points(df1$"P(+)"[df1$Design == design],df1$"mean SS"[df1$Design == design],type = "b", col = df1$col[df1$Design == design], pch = 16, cex = 2, lty = 2)
    }
    if(tested){
      points(df1$"P(+)"[df1$Design == design],df1$"mean tested"[df1$Design == design],type = "b", col = df1$col[df1$Design == design], pch = 16, cex = 2, lty = 2)
    }
  }

  axis(1, padj = - 0.3)
  axis(2, las=2, hadj = 1, padj = 0.25)
  if(xlab){
    mtext(expression(P(b^{"+"})), 1, line = 2.5, cex = 1.2)
  }
  if(ylab){
    if(plotMean){
      mtext("Mean sample size", 2, line = 3, cex = 1.5)
    } else {
      mtext("Mean screened", 2, line = 3, cex = 1.5)
    }
  }
  box()

}

plotTime = function(df1, ylab = TRUE, xlab = TRUE, ylim){
  par(mar=c(3,3,0.1,0.1))
  # plot(c(df1$"P(+)",df1$"P(+)"),c(df1[,"mean Time"],df1[,"95% Time"]),col=0, axes = FALSE, xlab = "", ylab = "")
  plot(df1$"P(+)",df1[,"mean Time"],col=0, axes = FALSE, xlab = "", ylab = "", ylim = ylim)
  for(design in the_design_names){
    points(df1$"P(+)"[df1$Design == design],df1[df1$Design == design,"mean Time"],type = "b", col = df1$col[df1$Design == design], pch = 16, cex = 2, lty = 2)

    #points(df1$"P(+)"[df1$Design == design],df1[df1$Design == design,"95% Time"],type = "b", col = df1$col[df1$Design == design], pch = 16, cex = 1.7, lty = 2)
  }

  axis(1, padj = - 0.3)
  axis(2, las=2, hadj = 1, padj = 0.25)
  if(xlab){
    mtext(expression(P(b^{"+"})), 1, line = 3.5, cex = 1.2)
  }
  if(ylab){
    mtext("Mean Trial time", 2, line = 3.7, cex = 1.5)
  }
  box()

}

plotTruePositive = function(df1, which, ylab = TRUE, xlab = FALSE, ylim = c(0.7,1), lineIfEqual = FALSE){

  par(mar=c(3,3,0.1,0.1))
  plot(df1[,"P(+)"],df1[, which],col=0, axes = FALSE, xlab = "", ylab = "", ylim = ylim)
  for(design in the_design_names){
    if(lineIfEqual && length(table(df1[df1$Design == design, which])) == 1) {
      abline(h=df1[df1$Design == design, which], col = df1$col[df1$Design == design], lwd = 2)
    } else {
    points(df1[df1$Design == design,"P(+)"],df1[df1$Design == design, which],type = "b", col = df1$col[df1$Design == design], pch = 16, cex = 2, lty = 2)
    }
  }

  axis(1, padj = - 0.3)
  axis(2, las=2, hadj = 1, padj = 0.25)
  if(xlab){
    mtext(expression(P(b^{"+"})), 1, line = 3.5, cex = 1.5)
  }
  if(ylab){
    mtext("Probability of correct conclusion", 2, line = 3.7, cex = 1.5)
  }
  box()

}
