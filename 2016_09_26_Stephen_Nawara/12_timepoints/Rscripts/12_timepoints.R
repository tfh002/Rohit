### Manage directories
setwd("D:/Tang/Rohit/2016_09_26_Stephen_Nawara/12_timepoints")
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
#goLookUp = oboToList(goLookUp)

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

# Calculate the slopes of divergence 
dat$DivergenceSlopes = NA
Divergence = t(abs(apply(dat[,tExpColmns], 1, normalize) - 
                     apply(dat[,pExpColmns], 1, normalize)))
NotNA = !is.na(Divergence[,1])
dat$DivergenceSlopes[NotNA] = apply(Divergence[NotNA,], 1, 
                                     function(x) coef(lm(x~Hrs))[2]) 


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
###### This gets us the 674 genes for which we have all data. 
geneNames = dat$Feature[!is.na(dat[,tExpColmns[1]]) & 
                          !is.na(dat[,pExpColmns[1]]) &
                          !is.na(dat$RLSMean)]

# Sort by RLS mean
dat = dat[order(dat$RLSMean, decreasing = T),]

# Putting only 674 genes for which we have everything from geneNames into dat
dat = dat[dat$Feature %in% geneNames,]



# Plot panels of selected genes
plotChoices(geneNames,  nr = 3, nc = 3, 
            makePDF = T, PDFdim = c(12, 9), 
            PDFname = "Levels_single.pdf",  mult = 1, 
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
                 PDFname = "Levels_multi.pdf",  mult = 1, 
                 lwd = 2,                  cex = 2, 
                 cex.main = 2, mar = c(5.1, 4.5, 6.0, 5.1), 
                 cex.x.axis = 0.8,         cex.y.axis = 2, 
                 pch1  = 24,               pch2 = 23,  
                 ylab1 = "mRNA Levels",    col1 = rainbow(length(geneNames)), 
                 ylab2 = "Protein Levels", col2 = "green", 
                 plotOnly = "tran", 
                 main = "Choices", addToPrev = F, norm = T, overTitle = "Overall Title")

plotComplexColmn("Pathways", "; ", "Pathways674.pdf")

#plotComplexColmn("GoTerm", "; ", "GoTerm.pdf")

# Finding the quantiles

