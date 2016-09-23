#setwd("D:/Tang/Rohit")

LS = read.csv("InputData/RLS.tsv", sep = "\t", stringsAsFactors = F, header = T)
LS = LS[!LS$Name=="",]
duplicatesMatt = LS$Name[duplicated(LS$Name)]
LSduplicates = LS[LS$Name %in%  duplicatesMatt,]

write.csv(LSduplicates, "duplicated gene names in Matts RLS file.csv", row.names = F)
LS2 = LS[!(LS$Name %in%  duplicatesMatt),]
dat = combine(dat,LS2,"Feature", "Name")
# add average column to the dat table

avgRLS = matrix(nrow = nrow(dat), ncol = 1)
for(i in 1:nrow(dat)){
  avgRLS[i] = mean(as.numeric(strsplit(dat$RLS[i], ",")[[1]]))
}
dat$buds = avgRLS

dat = dat[order(dat$buds, decreasing = T),]
rownames(dat) = 1:nrow(dat)

dat[1:100,c("Feature", "buds")]

#fix the NA's in dat

MattRLS = dat[,c("Feature", "buds", "RLS")]
head(MattRLS)
dim(MattRLS)
write.csv(MattRLS, "MattRLS3.csv", row.names = F)

######################################################################

#For each 5% fraction of the divergence histogram, how many of them are anti-longevity genes
#Criteria for anti-longevity genes
#Top 5%
#Top 10%
#Top 20%
#Top 30%




getwd()

dat[1:20,]
dim(LS)
head(dat)
dat_longevity = dat
write.csv(dat_longevity,"dat_longevity.csv", row.names = F)
getwd()
