setwd("D:/Tang/Rohit")
pdf("30 breaks DivergenceSlopeHistogram.pdf", width = 45, height = 45)
par(mar=c(10.1,8.1,8.1,2.1))
H = hist(dat2$Slopes, breaks = 30, col="yellow", xlab = "", ylab = "", main = "", xaxt = "n", yaxt = "n", ylim = c(0,120))
text(x=H$mids, y=H$counts, labels=paste(H$counts, " /", sum(H$counts), "=", round(H$counts/sum(H$counts),2)), pos = 3, cex = 2, col = "red")

axis(side = 1, at = H$mids, labels = H$mids, cex.axis = 2, line = -6, font = 2)
axis(side = 2, at = axTicks(2), labels = axTicks(2),cex.axis = 4, font = 10)
mtext("Slope of Divergence per Hour", side = 1,cex = 3)
mtext("Number of Genes", side = 2, cex = 3)
title(c("Linear Trend in Divergence between mRNA and Protein Expression Levels"), cex.main = 5, font = 2)



Quants = quantile(dat2$RLSMean, probs = 1 - c(0.05, 0.1, 0.2, 0.3))
for (i in 1:(length(H$breaks)-1)){
  InBin = H$breaks[i]<=dat2$Slopes & H$breaks[i+1]>dat2$Slopes
  lbls = dat2$Feature[InBin]
  RLS = dat2$RLSMean[InBin]
  Proport = round(sapply(Quants, function(x) length(which(RLS>x))/length(RLS)),2)
  text(H$mids[i], H$counts[i]*4/5, pos = 1, offset = 0, col="blue", font = 2, cex = 1.5, labels = paste(lbls, collapse = "\n"))
  text(H$mids[i],  H$counts[i], pos = 3, offset = 10, col= "darkgreen", cex = 2, font = 2, labels = paste(paste(names(Proport), Proport, sep = ": ")
                                                                                                          , collapse ="\n" ))
  
}
dev.off()

