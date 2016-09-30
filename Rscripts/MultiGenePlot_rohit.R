library(gridExtra)
library(ggthemes)

ScriptDir = "Rscripts"

source(paste0(ScriptDir,"/Main Slope Correlation.R"))


# Creating the Age Filter
Age_filter <- unique(Dat_pt$Age)[order(unique(Dat_pt$Age))%in%seq(3,10,1)]

# Creating the Gene Filter ranked by Average Life
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

# Creating the dataset needed to plot
Dat_plot <- Dat_pt%>%
  filter(Feature %in% Gene_filtered$Feature)%>%
  inner_join(y = GeneRank, by = c("Feature" = "Feature"))%>%
  select(c(gRowNum,Feature, ORF, Age,ProtExp,TranExp,ProtTranRatio))


# Generate a list of plots for each gene
p<-list()

for( i in GeneRank$gRowNum)
{
  Dat_s <- Dat_plot[Dat_plot$gRowNum == i , ]
  
  p[[i]]<- ggplot(Dat_s) + 
    geom_line(aes(x= Age, y=ProtTranRatio), col="green",size =2) +
    geom_point(aes(x= Age, y=ProtTranRatio),color = "red",fill = "blue", shape = 24, size =4) +
    geom_hline(aes(yintercept=1), colour="red", linetype=1, size =2) +
    geom_vline(aes(xintercept=Age_filter[1]), colour="blue", linetype="dashed", size =1) +
    geom_vline(aes(xintercept=Age_filter[8]), colour="blue", linetype="dashed", size =1)+
    ggtitle(paste(Dat_s$Feature[1],":",Dat_s$ORF[1])) + 
    ylab("Proteomic Transcriptomic Ratio") + 
    theme_gdocs() + 
    theme(axis.title.x = element_text(face="bold", color="black",size=12)) + 
    theme(axis.title.y = element_text(face="bold", color="black",size=12))
}



# Print the plots to a PDF file arranging them in 2x2 plots per page

pdf(paste0(plotDir, "/", "GenePlotsNorm_Rohit.pdf")) 

for( i in seq(0,length(p)-1,4))
{ grid.arrange(p[[i+1]],p[[i+2]],p[[i+3]],p[[i+4]], nrow =2, ncol=2, newpage = TRUE ) }

dev.off()



