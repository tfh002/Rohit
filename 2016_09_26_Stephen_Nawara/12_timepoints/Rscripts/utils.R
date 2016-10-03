### Helper Functions ###
clean   <- function(x){ 
    imp2Prime      = (x == "IMP2'")
    res            = tolower(gsub("([^[:alnum:] ])", "", x))
    res[imp2Prime] = "imp2'"
    return(res)
}

combine <- function(df1, df2, idColmn1, idColmn2){

	comp1 = clean(df1[,idColmn1])
	comp2 = clean(df2[,idColmn2])

	if(any(duplicated(comp1))){ print("Warning: Duplicates in df1") }
	if(any(duplicated(comp2))){ print("Warning: Duplicates in df2") }

    if(nrow(df1) > nrow(df2)){
        skp           = colnames(df2) == idColmn2
        res           = cbind(df1, df2[match(comp1, comp2), !skp, drop = F])
        extras        = df2[!(comp2 %in% clean(res[, idColmn1])), ]
        tmp           = data.frame(matrix(nrow = nrow(extras), ncol = ncol(res)))
        colnames(tmp) = colnames(res)

        colnames(extras)[colnames(extras) == idColmn2] = idColmn1
        tmp[, match(colnames(extras), colnames(res))]  = extras

        res = rbind(res, tmp)
    }else{
        skp = colnames(df1) == idColmn1
        res = cbind(df2, df1[match(comp2, comp1), !skp])
        extras        = df1[!(comp1 %in% clean(res[, idColmn2])), ]
        tmp           = data.frame(matrix(nrow = nrow(extras), ncol = ncol(res)))
        colnames(tmp) = colnames(res)

        colnames(extras)[colnames(extras) == idColmn1] = idColmn2
        tmp[, match(colnames(extras), colnames(res))]  = extras

        res = rbind(res, tmp)
    }

    rownames(res) = 1:nrow(res)
    return(res)
}

getExpColmns <- function(dat){
    tStart = grep("t7.8",   colnames(dat))
    tEnd   = grep("t72.3",  colnames(dat))
    pStart = grep("p7.8",   colnames(dat))
    pEnd   = grep("p72.3",  colnames(dat))

    tExpColmns = tStart:tEnd
    pExpColmns = pStart:pEnd
    return(list(tExpColmns = tExpColmns, 
                pExpColmns = pExpColmns, 
                Hrs        = as.numeric(gsub("t", "", 
                                             colnames(dat)[tExpColmns]))))
}



plotLevels <- function(Hrs, level, lwd = 2, pch = 24, cex = 2, mult = 1, 
                       cex.main = 2, mar = c(5.1, 4.5, 6.0, 5.1), 
                       cex.x.axis = 0.8, cex.y.axis = 2, ylab = "Levels", 
                       y.axis.side = 2, col = "blue", add = FALSE, main = "", 
                       panel = grid(col = "lightgray", lwd = 0.1, nx = NA, ny = NULL), 
                       addToPrev = F){
    #cex        = cex*mult
    #cex.main   = cex.main*mult
    #mar        = mar*mult
    lwd        = lwd*mult
    #cex.x.axis = cex.x.axis*mult
    #cex.y.axis = cex.y.axis*mult

    par(mar=mar)

    if(add == TRUE){ par(new=TRUE); panel = NULL; main = "" }
    plot(Hrs, level, type = "l", lwd = lwd, col = col, 
         yaxt = "n", xaxt = "n", ylab = "", xlab = "Hours",
         cex.lab = 1,  panel.first = panel)
    abline(v = Hrs, col = "lightgray", lwd = 0.1, lty = 2)
    points(Hrs, level, pch = pch, cex = cex, col = "red", bg = col)

    if(!addToPrev){
        axis(side = y.axis.side, cex.axis = cex.y.axis, col.axis = col, font=2)
        axis(side = 1, cex.axis = cex.x.axis, col.axis = "red", font = 2, at = Hrs, labels = Hrs)
        mtext(ylab, side = y.axis.side, line = 3, col = col, font = 2, cex = 1.3)

        for(i in 1:length(main)){
            title(main=main[i], line=5-1.5*i, cex.main=cex.main, font=2)
        }
    }
}



normalize <-function(x){ 
    if(is.null(nrow(x))){
        (x-min(x))/(max(x)-min(x)) 
    }else{
        t(apply(x, 1, function(r) (r-min(r))/(max(r)-min(r)) ))
    }
}






plotGene <- function(idx, lwd = 2, cex = 2, mult = 1,
                     cex.main = 2, mar = c(5.1, 4.5, 6.0, 5.1), 
                     cex.x.axis = 0.8,         cex.y.axis = 2, 
                     pch1  = 24,               pch2 = 23,  
                     ylab1 = "mRNA Levels",    col1 = "blue", 
                     ylab2 = "Protein Levels", col2 = "green", 
                     plotOnly = c("tran", "prot", "both")[3], 
                     main = "", addToPrev = FALSE, norm = FALSE, 
                     plotDiff = T){
    if(main == ""){
        main=c(paste0(idx, "  -  ", dat$Feature[idx], "   RLS: ", round(dat$RLSMean[idx],2)),  
               dat$Feature.Systematic.Name[idx], dat$Feature.Type[idx])
    }

    tranLevels = dat[idx, tExpColmns]
    protLevels = dat[idx, pExpColmns]

    if(plotDiff){
        diffLevels  = abs(normalize(tranLevels) - normalize(protLevels))
    }

    if(norm){ 
        tranLevels = normalize(tranLevels)
        protLevels = normalize(protLevels)
    }


    plotTran = ifelse(any(!is.na(tranLevels)) & plotOnly != "prot", TRUE, FALSE)
    plotProt = ifelse(any(!is.na(protLevels)) & plotOnly != "tran", TRUE, FALSE)
    errMsg   = paste0("\n", main[1], "- Nothing to plot\n")


    plot1 <- function(add = FALSE){
        plotLevels(Hrs, tranLevels, lwd = lwd, cex = cex, mult = mult, 
                   cex.main = cex.main, mar = mar, 
                   cex.x.axis = cex.x.axis, cex.y.axis = cex.y.axis, 
                   pch = pch1, ylab = ylab1, col = col1, 
                   main = main, add = add, addToPrev = addToPrev)
    }

    plot2 <- function(add = FALSE){
        plotLevels(Hrs, protLevels, lwd = lwd, cex = cex, mult = mult, 
                   cex.main = cex.main, mar = mar, 
                   cex.x.axis = cex.x.axis, cex.y.axis = cex.y.axis, 
                   pch = pch2, ylab = ylab2, col = col2, y.axis.side = 4, 
                   main = main, add = add, addToPrev = addToPrev)
    }

    plot3 <- function(add = FALSE){
        plotLevels(Hrs, diffLevels, lwd = lwd, cex = cex, mult = mult, 
                   cex.main = cex.main, mar = mar, 
                   cex.x.axis = cex.x.axis, cex.y.axis = cex.y.axis, 
                   pch = ".", ylab = "", col = "Red", y.axis.side = 4, 
                   main = main, add = add, addToPrev = T)
        abline(lm(as.numeric(diffLevels)~Hrs), col = "Orange")
    }

    if(!plotTran & !plotProt){ cat(errMsg) }

    if(!addToPrev){
        if(plotTran  & plotProt) { plot1(); plot2(TRUE); plot3(TRUE) }
        if(plotTran  & !plotProt){ plot1() }
        if(!plotTran & plotProt) { plot2() }
    }else{
        if(plotTran  & plotProt) { plot1(TRUE); plot2(TRUE); plot3(TRUE) }
        if(plotTran  & !plotProt){ plot1(TRUE) }
        if(!plotTran & plotProt) { plot2(TRUE) }

    }
}




plotChoices <- function(geneNames,  nr = 3, nc = 3, 
                        makePDF = F, PDFdim = c(12, 9),
                        PDFname = "levels.pdf",  mult = 1, 
                        lwd = 2,                  cex = 2, 
                        cex.main = 2, mar = c(5.1, 4.5, 6.0, 5.1), 
                        cex.x.axis = 0.8,         cex.y.axis = 2, 
                        pch1  = 24,               pch2 = 23,  
                        ylab1 = "mRNA Levels",    col1 = "blue", 
                        ylab2 = "Protein Levels", col2 = "green", 
                        plotOnly = c("tran", "prot", "both")[3], 
                        main = "", addToPrev = FALSE, norm = FALSE, 
                        overTitle = ""){

    PDFdim    = PDFdim*mult

    commonIdx = match(clean(geneNames), clean(dat$Feature))
    systemIdx = match(clean(geneNames), clean(dat$Feature.Systematic.Name))
    choiceIdx = unique(c(commonIdx, systemIdx))
    choiceIdx = choiceIdx[!is.na(choiceIdx)]

    if(makePDF){ 
        pdf(paste0(plotDir, "/", PDFname), width = PDFdim[1], height = PDFdim[2]) 
    }
    par(mfrow = c(nr, nc), oma=c(0,0,2,0))
    for(i in choiceIdx){
        plotGene(i, lwd = lwd, cex = cex, mult = mult, 
                 cex.main = cex.main, mar = mar, 
                 cex.x.axis = cex.x.axis, cex.y.axis = cex.y.axis, 
                 pch1  = pch1,  pch2 = pch2,  
                 ylab1 = ylab1, col1 = col1, 
                 ylab2 = ylab2, col2 = col2, 
                 plotOnly = plotOnly, 
                 main = main, addToPrev = addToPrev, norm = norm)
    }
    title(overTitle, outer = T, cex.main=2.5, col.main= "purple" )
    if(makePDF){ dev.off() }
}



plotChoicesMulti <- function(geneNames,  nr = 3, nc = 3, 
                             makePDF = F, PDFdim = c(12, 9),
                             PDFname = "levels.pdf",  mult = 1,
                             lwd = 2,                  cex = 2, 
                             cex.main = 2, mar = c(5.1, 4.5, 6.0, 5.1), 
                             cex.x.axis = 0.8,         cex.y.axis = 2, 
                             pch1  = 24,               pch2 = 23,  
                             ylab1 = "mRNA Levels",    col1 = "blue", 
                             ylab2 = "Protein Levels", col2 = "green", 
                             plotOnly = c("tran", "prot", "both")[3], 
                             main = "", addToPrev = FALSE, norm = FALSE, 
                             overTitle = ""){

    PDFdim    = PDFdim*mult

    commonIdx = match(clean(geneNames), clean(dat$Feature))
    systemIdx = match(clean(geneNames), clean(dat$Feature.Systematic.Name))
    choiceIdx = unique(c(commonIdx, systemIdx))
    choiceIdx = choiceIdx[!is.na(choiceIdx)]

    if(makePDF){ 
        pdf(paste0(plotDir, "/", PDFname), width = PDFdim[1], height = PDFdim[2]) 
    }
    #par(mfrow = c(nr, nc))
    par(oma=c(0,0,2,0))
    layout(matrix(c(1, 1, 1, 2), nrow = 4))
    plotGene(choiceIdx[1], lwd = lwd, cex = cex, mult = mult, 
             cex.main = cex.main, mar = mar, 
             cex.x.axis = cex.x.axis, cex.y.axis = cex.y.axis, 
             pch1  = pch1,  pch2 = pch2,  
             ylab1 = ylab1, col1 = col1[1], 
             ylab2 = ylab2, col2 = col2[1], 
             plotOnly = plotOnly, 
             main = main, addToPrev = F, norm = norm)
    for(i in choiceIdx[-1]){
        plotGene(i, lwd = lwd, cex = cex, mult = mult, 
                 cex.main = cex.main, mar = mar, 
                 cex.x.axis = cex.x.axis, cex.y.axis = cex.y.axis, 
                 pch1  = pch1,  pch2 = pch2,  
                 ylab1 = ylab1, col1 = col1[i], 
                 ylab2 = ylab2, col2 = col2[i], 
                 plotOnly = plotOnly, 
                 main = main, addToPrev = T, norm = norm)
    }
    title(overTitle, outer = T)
    par(mar=c(2.1, 4.5, 0.5, 5.1))
    plot(0, xaxt="n", yaxt="n", bty="n", pch="", ylab="", xlab="")
    legend("top", legend = dat$Feature[choiceIdx], col = col1, lwd = lwd*1.5, 
           cex = cex.x.axis*1.5, bty = "n", ncol = floor(length(choiceIdx)/2))

    if(makePDF){ dev.off() }
}


oboToList <- function(obo){

    starts = grep("[Term]", obo, fixed = T)
    ends   = c(starts[-1]-1, length(obo))

    res = list()
    for(i in 1:length(starts)){
        sub = obo[(starts[i] + 1):ends[i]]
        nms = unname(sapply(sub, function(x) unlist(strsplit(x, ":"))[1]))
        val = t(matrix(unname(sapply(sub, function(x) 
                                     paste(unlist(strsplit(x, ":"))[-1], collapse = ":")))))
        colnames(val) = nms
        goID          = unlist(strsplit(val[1], "GO:"))[2]
        res[[goID]] = as.data.frame(val)
        if(i %% 1000 ==  0){ print(paste(i, length(starts), sep = "/")) }
    }
    return(res)
}


plotComplexColmn <- function(colmn,  pattern,  PDFname = "Pathways.pdf"){

    uniqVals = unique(unlist(strsplit(dat[, colmn], pattern)))
    uniqVals = uniqVals[!is.na(uniqVals)]

    pdf(paste0(plotDir, "/", PDFname), width = 12,  height = 9)
    for(i in 1:length(uniqVals)){
        geneNames = dat$Feature[grep(uniqVals[i], dat[, colmn])]

        plotChoices(geneNames,  nr = 3, nc = 3, 
                    makePDF = F, PDFdim = c(12, 9), 
                    PDFname = "levels.pdf",  mult = 1, 
                    lwd = 2,                  cex = 2, 
                    cex.main = 2, mar = c(5.1, 4.5, 6.0, 5.1), 
                    cex.x.axis = 0.8,         cex.y.axis = 2, 
                    pch1  = 24,               pch2 = 23,  
                    ylab1 = "mRNA Levels",    col1 = "blue", 
                    ylab2 = "Protein Levels", col2 = "green", 
                    plotOnly = "both", 
                    main = "", addToPrev = F, norm = F, 
                    overTitle = uniqVals[i] )
        if(i %% 10 == 0){ print(i) }
    }
    dev.off()
}

