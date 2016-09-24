#setwd("D:/Tang/Rohit")

LS = read.csv("InputData/RLS.tsv", sep = "\t", stringsAsFactors = F, header = T)
LS = LS[!LS$Name=="",]
duplicatesMatt = LS$Name[duplicated(LS$Name)]
LSduplicates = LS[LS$Name %in%  duplicatesMatt,]

write.csv(LSduplicates, "duplicated gene names in Matts RLS file.csv", row.names = F)
LS2 = LS[!(LS$Name %in%  duplicatesMatt),]
LS3 = as.data.frame(sapply(LS2,toupper))

#dat = read.csv("dat file.csv", stringsAsFactors = F, header = T)
dat_inner = read.csv("dat inner file.csv", stringsAsFactors = F, header = T)



avgRLS = matrix(nrow = nrow(LS2), ncol = 1)
for(i in 1:nrow(LS2)){
  avgRLS[i] = mean(as.numeric(strsplit(LS2$RLS[i], ",")[[1]]))
}


LS3$buds = avgRLS

dat = merge(x = dat_inner,y = LS3,by.x = "Feature",by.y = "Name")



#dat = combine(dat,LS2,"Feature", "Name")
# add average column to the dat table

dat = dat[order(dat$buds, decreasing = T),]
rownames(dat) = 1:nrow(dat)


#########################################

MattRLS = dat[,c("Feature", "buds", "RLS","Sum","dev")]
MattRLS = MattRLS[order(MattRLS$dev,decreasing = F),]
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
