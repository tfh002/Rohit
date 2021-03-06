
# Loading Libraries
library("dplyr")
library("tidyr")
library("readr")
library("ggplot2")

# Setting up Directories
inputDir = "InputData"
plotDir  = "OutputPlots"        
fileDir  = "OutputFiles"
outDir   = "OutputData"

# Reading data from all sources
geneNames = read_csv(paste0(inputDir, "/InTransANDinProt_GeneNamesOnly.csv"))
prot      = read_csv(paste0(inputDir, "/Proteomic.csv"))
tran      = read_csv(paste0(inputDir, "/Transcriptomic.csv"))
gene      = read_csv(paste0(inputDir, "/AllChr.csv"))
path      = read_delim(paste0(inputDir, "/BiochemicalPathways.tab"),delim = "\t")
life      = read_tsv(paste0(inputDir, "/lifespan.tsv"))
goYeast   = read_delim(paste0(inputDir, "/gene_association.sgd"),delim = "\t")
RLS_df    = read_tsv(paste0(inputDir, "/RLS.tsv"),col_types = 'cc' )

# Cleaning data frames

prot_clean <- prot%>%
              gather(key = "Age", value = "ProtExp", - ORF)%>%
              mutate(Age = as.numeric(gsub("H","",Age)))

tran_clean <- tran%>%
              gather(key = "Age", value = "TranExp", - ORF)%>%
              mutate(Age = as.numeric(gsub("H","",Age)))

RLS_df <-  RLS_df[!RLS_df$Name=="",]

LS_clean  <- RLS_df %>%
             mutate(Name = toupper(Name))%>%
             separate(col = RLS,sep = ",",into = paste("V",c(1:250)),fill = "right")%>%
             gather(key = "colName",value = "Age", - Name,na.rm = TRUE)%>%
             select(-colName)%>%
             arrange(Name,Age)%>%
             mutate(Age = as.numeric(Age))%>%
             filter(!is.na(Age))%>%
             group_by(Name)%>%
             summarize(AvgLife = mean(Age))%>%
             arrange(desc(AvgLife))

life_clean <- life%>%
              select(Name = `Gene Symbol`, Longetivity = `Longevity Influence`)%>%
              filter(!is.na(Longetivity))%>%
              mutate(Anti = ifelse(Longetivity=="anti",1,0))%>%
              mutate(Pro = ifelse(Longetivity=="pro",1,0))%>%
              mutate(Fitness = ifelse(Longetivity=="fitness",1,0))%>%
              group_by(Name)%>%
              mutate(AntiProb = Anti/sum(Anti,Pro,Fitness))%>%
              summarize(AntiProb = round(sum(AntiProb),2))

# Joining data frames and calculating Prot Tran Ratio and correlations with Age
Dat_pt = inner_join(x = prot_clean,y = tran_clean, by = c("ORF" = "ORF", "Age" = "Age"))%>%
         mutate(ProtTranRatio = ProtExp/TranExp)%>%
         inner_join(y = gene,by = c("ORF" = "Feature.Systematic.Name"))%>%
         select(Feature,ORF,Age,ProtExp,TranExp,ProtTranRatio)%>%
         inner_join(y = LS_clean,by = c("Feature" = "Name"))%>%
         group_by(Feature)%>%
         mutate(Slope = cor(Age,ProtTranRatio))%>%
         ungroup(Feature)%>%
         mutate(GroupRank = ntile(Slope, 20))%>%
         arrange(Slope,ORF,Age)


Dat_pt_all_format = Dat_pt%>%
  select(ORF,Feature,AvgLife,Slope,Age,ProtTranRatio,GroupRank)%>%
  mutate(Age = paste0("H",Age))%>%
  spread(key = Age,value = ProtTranRatio)

# Preparing data frame for plotting histogram
Dat_hist = Dat_pt%>%
           inner_join(y = life_clean,by = c("Feature" = "Name"))%>%
           select(ORF,Feature,AvgLife,Slope,Age,ProtTranRatio,AntiProb)%>%
           mutate(Age = paste0("H",Age))%>%
           spread(key = Age,value = ProtTranRatio)

# Plotting the histogram
pdf(file = paste0(plotDir,"/Correlation Histogram.pdf"))

plot_hist <- ggplot(data = Dat_hist,mapping = aes(x = Slope))+
             geom_histogram(fill = "blue", color = "black",bins = 20)+
             labs(title = "Distribution of Age and Proteome/Transcriptome Correlation",
                  xlab = "Correlation")+
             facet_wrap(~ AntiProb)

print(plot_hist)
dev.off()

# Filtering histogram data for Anti Gene probability and cutting it into 5% groups
Data_Anti = Dat_hist%>%
            filter(AntiProb >= 0.5)%>%
            select(ORF,Feature,AvgLife,Slope,AntiProb)%>%
            ungroup(Feature)%>%
            mutate(GroupRank = ntile(Slope, 20))%>%
            arrange(GroupRank)

# Saving the data files
write.csv(Data_Anti, paste0(outDir,"/2nd ed Anti Gene Rank file.csv"), row.names = F)
write.csv(Dat_hist, paste0(outDir, "/2nd ed Histogram Data file.csv"), row.names = F)
write.csv(Dat_pt, paste0(outDir, "/2nd ed Proteomic Transcriptomic Data file.csv"), row.names = F)
write.csv(Dat_pt_all_format, paste0(outDir, "/2nd ed PT Data file 691 genes.csv"), row.names = F)
head(Dat_pt_all_format)

          
Dat_pt[1:100,]
head(Dat_pt_all_format)

