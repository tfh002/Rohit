TPdat = dat
Tdat = dat
head(TPdat)
head(Tdat)
colnames(Tdat)
Tdat<-Tdat[,-c(3:11,(ncol(Tdat)-1))]
Tdat<-Tdat
colnames(Tdat)
dim(Tdat)
Tdat<-Tdat[!is.na(Tdat$t7.8),]
head(Tdat)
Tdat<-Tdat[!is.na(Tdat$buds),]
Tdat<-Tdat[,-c(15:26)]
typoe.of(Tdat)
head(Tdat)
dim(Tdat)
Tdat<-Tdat[!is.na(Tdat$t7.8 & Tdat$buds),]
dim(Tdat)
colnames(Tdat)
head(Tdat)
write.csv(Tdat,"Tdat2.csv", row.names = F)
TPdat = Tdat
TPdat = TPdat[!is.na(TPdat$p7.8),]
dim(TPdat)
write.csv(Tdat, "Tdat.csv", row.names = F)
write.csv(TPdat, "TPdat.csv", row.names = F)
getwd()

TopTP0.3dat = TPdat

TopPdat<-function(data,percent){
TopPdat<-data[1:round((nrow(data)*percent/100),0),]
}

Top30pTPdat<-TopPdat(TPdat,30)
Top20pTPdat<-TopPdat(TPdat,20)
Top10pTPdat<-TopPdat(TPdat,10)
Top5pTPdat <-TopPdat(TPdat,5)

Top30pTdat<-TopPdat(Tdat,30)
Top20pTdat<-TopPdat(Tdat,20)
Top10pTdat<-TopPdat(Tdat,10)
Top5pTdat <-TopPdat(Tdat,5)


