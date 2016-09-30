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

# Slightly modified functions by Stephen to generate the plots 
mod_plotChoices <- function(geneNames,  nr = 3, nc = 3, 
                        makePDF = F, PDFdim = c(12, 9),
                        PDFname = "levels.pdf",  
                        lwd = 2,                  cex = 2, 
                        cex.x.axis = 0.8,         cex.y.axis = 2, 
                        pch1  = 24,               pch2 = 23,  
                        ylab1 = "mRNA Levels",    col1 = "blue", 
                        ylab2 = "Protein Levels", col2 = "green", 
                        plotOnly = c("tran", "prot", "both")[3], 
                        main = "", addToPrev = FALSE, norm = FALSE){
  
  # removed clean up and matching code here to just directly include index of choices
  choiceIdx = as.integer(row.names(geneNames))
  
  if(makePDF){ 
    pdf(paste0(plotDir, "/", PDFname), width = PDFdim[1], height = PDFdim[2]) 
  }
  par(mfrow = c(nr, nc))
  for(i in choiceIdx){
    mod_plotGene(i, lwd = lwd, cex = cex, 
             cex.x.axis = cex.x.axis, cex.y.axis = cex.y.axis, 
             pch1  = pch1,  pch2 = pch2,  
             ylab1 = ylab1, col1 = col1, 
             ylab2 = ylab2, col2 = col2, 
             plotOnly = plotOnly, 
             main = main, addToPrev = addToPrev, norm = norm)
  }
  if(makePDF){ dev.off() }
}


mod_plotGene <- function(idx, lwd = 2, cex = 2, 
                     cex.x.axis = 0.8,         cex.y.axis = 2, 
                     pch1  = 24,               pch2 = 23,  
                     ylab1 = "mRNA Levels",    col1 = "blue", 
                     ylab2 = "Protein Levels", col2 = "green", 
                     plotOnly = c("tran", "prot", "both")[3], 
                     main = "", addToPrev = FALSE, norm = FALSE){
  
  # Modified function to accept Dat_plot instead of dat which was being used
  # Used Dat_plot to get relevant data for single gene using only the idx
  
  Dat_single = Dat_plot[Dat_plot$gRowNum == idx,]
  
  if(main == ""){
    main=c(paste(idx, Dat_single$Feature[1], sep = ": "),Dat_single$ORF[1], "ORF")
  }
  
  tranLevels = Dat_single$TranExp
  protLevels = Dat_single$ProtExp
  Hrs        = Dat_single$Age
  
  if(norm){ 
    tranLevels = normalize(tranLevels)
    protLevels = normalize(protLevels)
  }
  
  plotTran = ifelse(any(!is.na(tranLevels)) & plotOnly != "prot", TRUE, FALSE)
  plotProt = ifelse(any(!is.na(protLevels)) & plotOnly != "tran", TRUE, FALSE)
  errMsg   = paste0("\n", main[1], "- Nothing to plot\n")
  
  
  plot1 <- function(add = FALSE){
    plotLevels(Hrs, tranLevels, lwd = lwd, cex = cex, 
               cex.x.axis = cex.x.axis, cex.y.axis = cex.y.axis, 
               pch = pch1, ylab = ylab1, col = col1, 
               main = main, add = add, addToPrev = addToPrev)
  }
  
  plot2 <- function(add = FALSE){
    plotLevels(Hrs, protLevels, lwd = lwd, cex = cex, 
               cex.x.axis = cex.x.axis, cex.y.axis = cex.y.axis, 
               pch = pch2, ylab = ylab2, col = col2, y.axis.side = 4, 
               main = main, add = add, addToPrev = addToPrev)
  }
  
  if(!plotTran & !plotProt){ cat(errMsg) }
  
  if(!addToPrev){
    if(plotTran  & plotProt) { plot1(); plot2(TRUE) }
    if(plotTran  & !plotProt){ plot1() }
    if(!plotTran & plotProt) { plot2() }
  }else{
    if(plotTran  & plotProt) { plot1(TRUE); plot2(TRUE) }
    if(plotTran  & !plotProt){ plot1(TRUE) }
    if(!plotTran & plotProt) { plot2(TRUE) }
    
  }
}


# Did not modify this function at all

plotLevels <- function(Hrs, level, lwd = 2, pch = 24, cex = 2, 
                       cex.x.axis = 0.8, cex.y.axis = 2, ylab = "Levels", 
                       y.axis.side = 2, col = "blue", add = FALSE, main = "", 
                       panel = grid(col = "lightgray", lwd = 0.1, nx = NA, ny = NULL), 
                       addToPrev = F){
  par(mar=c(5.1, 4.5, 6.0, 5.1))
  
  
  if(add == TRUE){ par(new=TRUE); panel = NULL; main = "" }
  plot(Hrs, level, type = "l", lwd = lwd, col = col, 
       yaxt = "n", xaxt = "n", ylab = "", xlab = "Hours",
       panel.first = panel)
  abline(v = Hrs, col = "lightgray", lwd = 0.1, lty = 2)
  points(Hrs, level, pch = pch, cex = cex, col = "red", bg = col)
  
  if(!addToPrev){
    axis(side = y.axis.side, cex.axis = cex.y.axis, col.axis = col, font=2)
    axis(side = 1, cex.axis = cex.x.axis, col.axis = "red", font = 2, at = Hrs, labels = Hrs)
    mtext(ylab, side = y.axis.side, line = 3, col = col, font = 2, cex = 1.3)
    
    for(i in 1:length(main)){
      title(main=main[i], line=5-1.5*i, cex.main=2, font=2)
    }
  }
}

# This is used in plot. 

normalize <-function(x){ (x-min(x))/(max(x)-min(x)) }