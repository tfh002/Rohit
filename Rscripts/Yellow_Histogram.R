setwd("D:/Tang/Rohit/2016_09_26_Stephen_Nawara/12_timepoints")
### Manage directories
#setwd("~/Thomas/12_timepoints")
inputDir = "InputData"
plotDir  = "OutputPlots" 
fileDir  = "OutputFiles"

require(digitize)
source("Rscripts/utils.R")

#geneNames = c(dat$Feature[1:4], dat$Feature.Systematic.Name[5:9])
geneNames<-read.csv(paste0(inputDir, "/InTransANDinProt_GeneNamesOnly.csv"), stringsAsFactors=F)[,]
### Read in data ###
# Main data
prot = read.csv(paste0(inputDir, "/Proteomic.csv"), stringsAsFactors=F)
tran = read.csv(paste0(inputDir, "/Transcriptomic.csv"), stringsAsFactors=F)

# Additional Info
gene     = read.csv(paste0(inputDir, "/AllChr.csv"),  stringsAsFactors = F)
path     = read.csv(paste0(inputDir, "/BiochemicalPathways.tab"),
                    stringsAsFactors = F, sep = "\t")
life     = read.csv(paste0(inputDir, "/lifespan.tsv"),
                    stringsAsFactors = F, sep = "\t")
goYeast  = read.csv(paste0(inputDir, "/gene_association.sgd"), sep = "\t",
                    stringsAsFactors = F)
goLookUp = scan(paste0(inputDir, "/go.obo"), what = "", sep = "\t")

RLS      =  read.csv(paste0(inputDir, "/yeast_RLS_data.txt"), sep = "\t",
                     stringsAsFactors = F,  header = F)



### Pre-processing ###
colnames(tran)[-1] = gsub("H", "t", colnames(tran)[-1])
colnames(prot)[-1] = gsub("H", "p", colnames(prot)[-1])

prot = prot[!grepl(";", prot$ORF),]

# Get proportion of references claiming pro, anti, fitness 
# for each gene in the lifespan data
lGenes           = unique(life$Gene.Symbol)
lProbs           = data.frame(matrix(nrow=length(lGenes), ncol=4))
colnames(lProbs) = c("CommonName", "pro", "anti", "fitness")

for(i in 1:length(lGenes)){
  sub            = life[life$Gene.Symbol == lGenes[i], ]
  Pr             = sapply(c("pro", "anti", "fitness"), function(x) 
    sum(sub$Longevity.Influence == x)/nrow(sub))
  lProbs[i, 1]   = lGenes[i]
  lProbs[i, 2:4] = Pr
}
lProbs[,2:4]     = sapply(lProbs[, 2:4], as.numeric)


# Arrange the pathway data so there is one row per gene
pathGenes = unique(clean(path$Gene))
uniqPaths = as.data.frame(matrix(nrow = length(pathGenes), ncol = 3))
for(i in 1:length(pathGenes)){
  sel             = path$Pathway[clean(path$Gene) == clean(pathGenes[i])]
  uniqPaths[i, 1] = pathGenes[i]
  uniqPaths[i, 2] = length(sel)
  uniqPaths[i, 3] = paste(sel, collapse = "; ")
}
colnames(uniqPaths) = c("Gene", "Npathways", "Pathways")


# Get GO IDs for each gene
geneGO = aggregate(GoTerm~Gene, goYeast, function(x) 
  gsub("GO:", "", paste(x, collapse = "; ")))

# Convert the GO obo file to R list
goLookUp = oboToList(goLookUp)

# Give RLS data colnames, remove duplicates (possibly combine duplicates later)
# Also extract means as a summary stat (maybe do more with this later)
colnames(RLS) = c("Name", "RLS")
RLS           = RLS[!duplicated(RLS$Name), ]

RLS$RLSMean = sapply(RLS$RLS, function(x) mean(as.numeric(strsplit(x,  ",")[[1]])))



# Combine the various data sources into one data frame
dat = combine(tran, prot,      "ORF",     "ORF")
dat = combine(dat,  gene,      "ORF",     "Feature.Systematic.Name")
dat = combine(dat,  uniqPaths, "Feature", "Gene")
dat = combine(dat,  lProbs,    "Feature", "CommonName")
dat = combine(dat,  geneGO,    "Feature", "Gene")
dat = combine(dat,  RLS,       "Feature", "Name")


# Remove undesirable rows/columns
# Determine the columns containing expression data
expColmns  = getExpColmns(dat)
tExpColmns = expColmns$tExpColmns
pExpColmns = expColmns$pExpColmns
Hrs        = expColmns$Hrs


# Remove rows without any expression data
remRows = apply(dat[, c(tExpColmns, pExpColmns)], 1, function(x) all(is.na(x)))
dat     = dat[!remRows, ]


# Reorder columns so that expression data is last
dat = dat[, c(1:(tExpColmns[1] - 1), 
              (pExpColmns[length(Hrs)] + 1):ncol(dat), 
              tExpColmns, pExpColmns) ]

# Remove columns not being used (for now)
remCols = sapply(dat, function(x) 
  length(unique(x)))  == 1 | colnames(dat) %in% c("Npathways")
dat     = dat[, !remCols]

# Columns have been reordered so get expression columns again
expColmns  = getExpColmns(dat)
tExpColmns = expColmns$tExpColmns
pExpColmns = expColmns$pExpColmns



## Make the plots



# Plot panels of selected genes
plotChoices(geneNames,  nr = 3, nc = 3, 
            makePDF = T, PDFdim = c(12, 9), 
            PDFname = "levels.pdf",  mult = 1, 
            lwd = 2,                  cex = 2, 
            cex.main = 2, mar = c(5.1, 4.5, 6.0, 5.1), 
            cex.x.axis = 0.8,         cex.y.axis = 2, 
            pch1  = 24,               pch2 = 23,  
            ylab1 = "mRNA Levels",    col1 = "blue", 
            ylab2 = "Protein Levels", col2 = "green", 
            plotOnly = "both", 
            main = "", addToPrev = F, norm = F,  overTitle = "Overal Title")

# Plot selected genes on one chart
plotChoicesMulti(geneNames,  nr = 3, nc = 3, 
                 makePDF = T, PDFdim = c(12, 9), 
                 PDFname = "levels2.pdf",  mult = 1, 
                 lwd = 2,                  cex = 2, 
                 cex.main = 2, mar = c(5.1, 4.5, 6.0, 5.1), 
                 cex.x.axis = 0.8,         cex.y.axis = 2, 
                 pch1  = 24,               pch2 = 23,  
                 ylab1 = "mRNA Levels",    col1 = rainbow(length(geneNames)), 
                 ylab2 = "Protein Levels", col2 = "green", 
                 plotOnly = "tran", 
                 main = "Choices", addToPrev = F, norm = T, overTitle = "Overall Title")

plotComplexColmn("Pathways", "; ", "Pathways.pdf")

#plotComplexColmn("GoTerm", "; ", "GoTerm.pdf")

# Analysis


dat2 = dat[!is.na(dat[,tExpColmns[1]]) & !is.na(dat[,pExpColmns[1]]) & !is.na(dat$RLSMean), ]
#dat2 = dat[!is.na(dat[,tExpColmns[1]]) & !is.na(dat[,pExpColmns[1]]),] 


# Calculate slope of the divergence: 

dat2$Slopes = sapply(1:nrow(dat2), function(i)
  coef(lm(as.numeric(abs(normalize(dat2[i,tExpColmns])-normalize(dat2[i,pExpColmns])))~Hrs))[2])




length(H$breaks)
dat$Feature[which(dat2$Slopes>0.015)]


pdf("DivergenceSlopeHistogram.pdf", width = 45, height = 45)
par(mar=c(10.1,8.1,8.1,2.1))
H = hist(dat2$Slopes, breaks = 20, col="yellow", xlab = "", ylab = "", main = "", xaxt = "n", yaxt = "n", ylim = c(0,120))
text(x=H$mids, y=H$counts, labels = H$counts, pos = 3, cex = 3, col = "red", font = 2)
axis(side = 1, at = H$mids, labels = H$mids, cex.axis = 5, line = 10)
axis(side = 2, at = axTicks(2), labels = axTicks(2))
mtext("Slope of Divergence per Hour", side = 1)
mtext("Number of Genes", side = 2)
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
  text(H$mids[i], H$counts[i], pos = 3, offset = 20, col = "orange", cex = 2, font = 2, labels = paste(length(H$counts)/674, collapse = "\n"))
  
}
dev.off()
write.csv(dat2, "dat2Hist.csv", row.names = F)

which(RLS<Quants[1])


#dat2[InBin, c("Feature", "Slopes")]


Sp = sapply(1:length(tExpColmns), function(i)
  cor(dat2[,tExpColmns[i]],dat2[,pExpColmns[i]], method = "spearman"))

plot(Hrs, Sp, ylim = c(0.6, 0.8))
abline(h=c(0.75, 0.7, 0.73), lty=2)

# Digitizing Janssens plot
if(FALSE){
  fName = paste0(inputDir, "/Janssens_uncoupling_correlation.png")
  cal=ReadAndCal(fName)
  pts = DigitData()
  JanssenSp = Calibrate(pts, cal, x1=6, x2 = 72, y1 = 0.65, y2 = 0.8)
  JanssenSp
  JanssenSp$x = Hrs
  JanssenSp
  
  Comp = cbind(JanssenSp, Sp)
  colnames(Comp) = c("Hrs", "JanssenSp", "Sp")
}


dat2$CorSP = sapply(1:nrow(dat2), function(i)
  cor(as.numeric(dat2[i,tExpColmns]), as.numeric(dat2[i,pExpColmns]), method = "spearman"))



dat3 = dat2[order(dat2$CorSP, decreasing = T),]
dim(dat3)
rownames(dat3)= 1:nrow(dat3)
head(dat3)
plot(as.numeric(dat3$RLSMean), as.numeric(dat3$CorSP))
head(dat3)
dat3[1,]
plot(dat3[1,tExpColmns], dat3[1, pExpColmns])
dat3[1:5,]
tail(dat3)
dat4 = dat3[,c("Feature", "RLSMean", "CorSP")]
colnames(dat3)
head(dat4)

plotGene(1)
dim(dat4)
dat3
plotChoices("UBP1", nr=1, nc=1)
tail(dat3)
which.min(abs(dat3$CorSP))
dat3$Feature[317]
CorSP_0 = dat4[order(abs(dat3$CorSP)),]

head(CorSP_0,100)
mean(dat3$RLSMean)
mean(dat3$RLSMean[1:100])
mean(rev(dat3$RLSMean)[1:100])
mean(CorSP_0$RLSMean[1:100])
tail(dat4,100)
head(dat4,100)
head(dat3)

head(dat3)
dat3$Slopes = sapply(1:nrow(dat3), function(i)
  coef(lm(as.numeric(abs(normalize(dat3[i,tExpColmns])-normalize(dat3[i,pExpColmns])))~Hrs))[2])
head(dat3)
dat3 = dat3[order(dat3$Slopes, decreasing = T),]
head(dat3)
require(rgl)
plot3d(dat3$RLSMean, dat3$CorSP, dat3$Slopes)
dat4
head(dat4)
dat4[dat4$RLSMean>30, ]
