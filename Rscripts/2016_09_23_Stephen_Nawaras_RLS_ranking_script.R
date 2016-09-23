setwd("D:/Tang/Rohit")

LS = read.csv("InputData/RLS.tsv", sep = "\t", stringsAsFactors = F, header = T)
LS = LS[!LS$Name=="",]
duplicatesMatt = LS$Name[duplicated(LS$Name)]
dupliRowsMatt =  LS[LS$Name[duplicated(LS$Name)],}
write.csv(dupliRowsMatt, "dupliRowsMatt.csv", row.names = F)
LS = LS[!duplicated(LS$Name),]
dat = combine(dat,LS,"Feature", "Name")

buds = matrix(nrow = nrow(dat), ncol = 1)
for(i in 1:nrow(dat)){
  buds[i] = mean(as.numeric(strsplit(dat$RLS[i], ",")[[1]]))
}
dat$buds = buds  

dat = dat[order(dat$buds, decreasing = T),]
rownames(dat) = 1:nrow(dat)

dat[1:1000,c("Feature", "buds")]

MattRLS = dat[,c("Feature", "buds", "RLS")]
head(MattRLS)
dim(MattRLS)
write.csv(MattRLS, "MattRLS3.csv", row.names = F)

######################################################################



getwd()

dat[1:20,]
dim(LS)
head(dat)
dat_longevity = dat
write.csv(dat_longevity,"dat_longevity.csv", row.names = F)
getwd()
