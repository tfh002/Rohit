########333
setwd()
data = dat
dim(data)
data1 = subset(data, select = -c(3:11))
dim(data1)
colnames(data1)
head(data1)
data2 = data1[complete.cases(data1),]
dim(data2)
head(data2)

tMin<-rep(NA,length(data2[,2])) # creates positions in which tMin can be put
tMax<-rep(NA,nrow(data2))
pMin<-rep(NA,nrow(data2))
pMax<-rep(NA,nrow(data2))

for (i in 1:nrow(data2)){
  tMin[i] = min(data2[i,(3:14)])
  tMax[i] = max(data2[i,(3:14)])
  pMin[i] = min(data2[i,(15:ncol(data2))])
  pMax[i] = max(data2[i,((ncol(data2)-11):ncol(data2))])
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



ord = order(data6[,"Sum"])
coupling = data6[ord,c("Feature", "Sum")]
rownames(coupling)=1:nrow(coupling)

geneNames = coupling$Feature
