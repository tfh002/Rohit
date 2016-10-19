########333
setwd()
data = TPdat
dim(data)
head(data)
data2 = data
head(data)
geneNames = data
tMin<-rep(NA,length(data2[,2])) # creates positions in which tMin can be put
tMax<-rep(NA,nrow(data2))
pMin<-rep(NA,nrow(data2))
pMax<-rep(NA,nrow(data2))

for (i in 1:nrow(data2)){
  tMin[i] = min(data2[i,(3:14)])
  tMax[i] = max(data2[i,(3:14)])
  pMin[i] = min(data2[i,(15:ncol(data2)-1)])
  pMax[i] = max(data2[i,((ncol(data2)-11):ncol(data2)-1)])
}# determines tMin for each line


colnames(data2)
head(data2)
data3<-data.frame(data2,tMin, tMax, tRange=tMax-tMin, pMin, pMax, pRange=pMax-pMin)
head(data3)

data4<-data.frame(matrix(NA,nrow=nrow(data3),ncol=ncol(data3)))

data4Colnames<-rep(NA,ncol(data4))

for(i in 1:ncol(data4)){
  if (i<3 | i>26) data4[,i]=data3[,i]
  if(i>=3 & i<=14) data4[,i]=data3[,i]-tMin
  if(i>=15 & i<=26) data4[,i]=data3[,i]-pMin
  
  if (i<3 | i>26) data4Colnames[i]=colnames(data3[i])
  if(i>=3 & i<=14)data4Colnames[i]=paste0(colnames(data3[i]),"_tMin")
  if(i>=15 & i<=26)data4Colnames[i]=paste0(colnames(data3[i]),"_pMin")
}
colnames(data4)<-data4Colnames
head(data4)




data5<-data.frame(matrix(NA,nrow=nrow(data4),ncol=ncol(data4)))

data5Colnames<-rep(NA,ncol(data4))

for(i in 1:ncol(data4)){
  if (i<3 | i>26) data5[,i]=data4[,i]
  if(i>=3 & i<=14) data5[,i]=data4[,i]/data4$tRange
  if(i>=15 & i<=26) data5[,i]=data4[,i]/data4$pRange
  
  if (i<3 | i>26) data5Colnames[i]=colnames(data3[i])
  if(i>=3 & i<=14)data5Colnames[i]=paste0(colnames(data3[i]),"n")
  if(i>=15 & i<=26)data5Colnames[i]=paste0(colnames(data3[i]),"n")
}
colnames(data5)<-data5Colnames
head(data5)


data6<-data.frame(matrix(NA,nrow=nrow(data5),ncol=ncol(data5)-12))
data6Colnames<-rep(NA,ncol(data4)-12)

for(i in 1:ncol(data6)){
  if (i<3) data6[,i]=data5[,i]
  if(i>=3 & i<=14) data6[,i]=data5[,i]-data5[,(i+12)]
  if(i>=15) data6[,i]=data5[,(i+12)]
  
  if (i<3) data6Colnames[i]=colnames(data3[i])
  if(i>=3 & i<=14)data6Colnames[i]=paste0(colnames(data3[i]), "d")
  if(i>=15)data6Colnames[i]=paste0(colnames(data3[i+12]),"d")
}

colnames(data6)<-data6Colnames
head(data6)
data6[,(3:20)]<-abs(data6[,(3:20)])
data4<-format(data4, scientific=FALSE)

Sum<-rowSums(data6[,3:14])
data6<-cbind(data6,Sum)
head(data6)

normalize(dat[,tExpColmns])-normalize(dat[,pExpColmns])
colnames(dat)[tExpColmns]
normalize(dat[2,tExpColmns])
nTrans = normalize(dat[,tExpColmns])
head(nTrans)
head(dat)



difSum = data.frame(Feature = dat$Feature)
tpDif = t(apply(dat[,tExpColmns], 1, normalize)) - t(apply(dat[,pExpColmns], 1, normalize))
difSum$Sum = rowSums(abs(tpDif))
difSum = difSum[order(difSum$Sum),]
difSum = difSum[!is.na(difSum$Sum),]



head(difSum)
head(dat)


head(data6)
ord = order(data6[,"Sum"])
coupling = data6[ord,c("Feature", "Sum", "budsd")]
rownames(coupling)=1:nrow(coupling)

geneNames = coupling$Feature
list675<-match(geneNames,data6$Feature)
dat675<-TPdat[list675,]
dat675<-dat675[order(dat675[,"buds"],decreasing = T),]


write.csv(dat675,"dat675.csv",row.names=F)


## You have to run plot function script for ploting rest output
dat=dat675

# Columns have been reordered so get expression columns again
expColmns  = getExpColmns(dat)
tExpColmns = expColmns$tExpColmns
pExpColmns = expColmns$pExpColmns

difSum = data.frame(Feature = dat$Feature, stringsAsFactors = F)
tpDif = t(apply(dat[,tExpColmns], 1, normalize)) - t(apply(dat[,pExpColmns], 1, normalize))


difSum$Sum = rowSums(abs(tpDif))
tpDif2 = tpDif[order(difSum$Sum),]
difSum = difSum[order(difSum$Sum),]
difSum = difSum[!is.na(difSum$Sum),]
tpDif2 = tpDif2[!is.na(tpDif2[,1]),]

dev = apply(tpDif2,1,function(x) coef(lm(abs(x)~Hrs))[2])
ord.dev = order(dev)
difSum = difSum[ord.dev,]
tpDif2 = tpDif2[ord.dev,]


pdf(paste0(plotDir, "/Trans-Prot_abs4_HistSahidul.pdf"), width = 12, height = 9)
hist(dev, breaks = 20, col="brown", xlab = "Slope: Change Per Hour", main = c("Linear Trend in Divergence", "Between mRNA and Protein Expression Levels"))
par(mfrow=c(3,3))
for(i in 1:nrow(tpDif2)){
  plotLevels(Hrs, abs(tpDif2[i,]), main = c(difSum$Feature[i],"Trans-Prot Magnitude"), col="red")
  #abline(h=0,lty=1,col="purple",lwd=3)
  abline(lm(abs(tpDif2[i,])~Hrs), lty=1, col= "orange", lwd=3)
}
dev.off()

pdf("2nd ed Divergence_versus_knockout_life_extension.pdf", width = 12, height = 9)
plot(dat2$Slopes, dat2$RLSMean, xlab = "Slope of Divergence per Hour", ylab = "Average RLS for knockout", main = "Rate of Divergence vs. Lifespan of Knockout")
abline(v=c(-.008, -0.006, 0.004, 0.006), h = Quants, col="red", lty=2)
text(rep(-0.015, 4), Quants, names(Quants), font = 2)
text(c(-0.007, 0.005), rep(5,2), "Anti-Life", font = 2)
dev.off()



ExpColumns = getExpColmns(dat2)
tExpColmns = ExpColumns$tExpColmns
TransDifs = t(apply(dat2[, tExpColmns], 1, diff)/diff(Hrs))
pExpColmns = ExpColumns$pExpColmns
ProtDifs = t(apply(dat2[, pExpColmns], 1, diff)/diff(Hrs))

ExpDifs = TransDifs - ProtDifs
DifSum = rowSums(ExpDifs)
absDifSum = rowSums(abs(ExpDifs))

DifsTransProt = dat2[, tExpColmns]-dat2[, pExpColmns]

colnames(DifsTransProt) = paste0(colnames(DifsTransProt), "Trans_minus_Prot_value_table")
colnames(TransDifs) = paste0(colnames(TransDifs), "TransDerivative")
colnames(ProtDifs) = paste0(colnames(ProtDifs), "ProtDerivative")
colnames(ExpDifs) = paste0(colnames(ExpDifs), "Trans_minus_Prot_Derivative")

file = cbind(dat2, TransDifs, ProtDifs, ExpDifs, DifSum, absDifSum, DifsTransProt)
head(file)

write.csv(file, "dataFileForExcel.csv", row.names = F)


dat2[1, tExpColmns]
dim(TransDifs)
dim(ProtDifs)
colnames(dat2)
getExpColmns(dat2)
getExpColmns(dat)
colnames(dat)
