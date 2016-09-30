
ScriptDir = "Rscripts"

source(paste0(ScriptDir,"/Utils.R"))
source(paste0(ScriptDir,"/Main Slope Correlation.R"))

# Querying DAT_PT for Age>=14 and when Proteomic Expression exceeds Transcriptomic Expression atleast once
Age_filter <- unique(Dat_pt$Age)[order(unique(Dat_pt$Age))%in%seq(3,10,1)]


Gene_filtered <- Dat_pt%>%
  filter(Age %in% Age_filter)%>%
  filter(ProtTranRatio>1)%>%
  mutate(Age = paste0("H",Age))%>%
  select(-c(ProtExp,TranExp))%>%
  spread(key = Age,value = ProtTranRatio)%>%   # spread to order the set in avglife
  arrange(desc(AvgLife))%>%
  select(Feature)

write.csv(Gene_filtered, paste0(outDir,"/Filtered Gene List.csv"), row.names = F)

GeneRank <- Gene_filtered%>%
  mutate(gRowNum = as.integer(row.names(Gene_filtered)))

Dat_plot <- Dat_pt%>%
  filter(Feature %in% Gene_filtered$Feature)%>%
  inner_join(y = GeneRank, by = c("Feature" = "Feature"))%>%
  select(c(gRowNum,Feature, ORF, Age,ProtExp,TranExp))
  

# Using Stephen's method of plotting with minor modifications

mod_plotChoices(Gene_filtered,  nr = 3, nc = 3, 
            makePDF = T, PDFdim = c(12, 9), 
            PDFname = "GenePlots.pdf",   
            lwd = 2,                cex = 2, 
            cex.x.axis = 0.8,         cex.y.axis = 2, 
            pch1  = 24,               pch2 = 23,  
            ylab1 = "mRNA Levels",    col1 = "blue", 
            ylab2 = "Protein Levels", col2 = "green", 
            plotOnly = "both", 
            main = "", addToPrev = F, norm = T)
