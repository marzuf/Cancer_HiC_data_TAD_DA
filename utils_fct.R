
printAndLog <- function(txt, logFile=""){
  cat(txt)
  cat(txt, file = logFile, append=T)
}

######################################################################################################################################################################################################
###################################################################################################################################################################################################### add_curv_fit (function)
######################################################################################################################################################################################################


densplot <- function(x,y, pch=19, cex=1, ...){
	df <- data.frame(x,y)
	d <- densCols(x,y, colramp=colorRampPalette(c("black", "white")))
	df$dens <- col2rgb(d)[1,] + 1L
	cols <- colorRampPalette(c("#000099", "#00FEFF", "#45FE4F","#FCFF00", "#FF9400", "#FF3100"))(256)
	df$col <- cols[df$dens]
	df <- df[order(df$dens),]
	plot(df$x,df$y, pch=pch, col=df$col, ...)
}


######################################################################################################################################################################################################
###################################################################################################################################################################################################### add_curv_fit (function)
######################################################################################################################################################################################################

add_curv_fit <- function(x, y, withR2 = TRUE, R2shiftX = 0, R2shiftY = 0,...) {
  mymodel <- lm(y~x)
  abline(mymodel, ...)
  if(withR2) {
    r2Txt <- paste0("adj. R2 = ", sprintf("%.2f", summary(mymodel)$adj.r.squared))
    r2X <- x[which.min(x)] + R2shiftX
    r2Y <- fitted(mymodel)[which.min(x)]
    text(x = r2X, y = r2Y, 
         labels = r2Txt, 
         adj=c(1,0),
         pos=3,
         cex = 0.7)
  }
}

######################################################################################################################################################################################################
###################################################################################################################################################################################################### addCorr (function)
######################################################################################################################################################################################################

addCorr <- function(x, y, legPos="topright", corMet="pearson", ...) {
  corMet <- tolower(corMet)
  stopifnot(corMet %in% c("pearson", "kendall", "spearman"))
  x_new <- x[!is.na(x) & !is.na(y)]
  y_new <- y[!is.na(x) & !is.na(y)]
  x <- x_new
  y <- y_new
  stopifnot(length(x) == length(y))

  if(length(x) < 3) {
    legTxt <- paste0(paste0(toupper(substr(corMet,1,1)), "CC"), " = NA", "\n", "(# obs. < 3)")
  } else {
    ct <- cor.test(x,y, method = corMet)
    corCoeff <- ct$estimate
    corPval <- ct$p.value
    legTxt <- paste0(paste0(toupper(substr(corMet,1,1)), "CC"), " = ", round(corCoeff, 4), "\n", "(p-val = ", sprintf("%2.2e", corPval), ")")
  }
  legend(legPos, legend = legTxt, ...)
}

######################################################################################################################################################################################################
###################################################################################################################################################################################################### plot_multiDens(function)
######################################################################################################################################################################################################

plot_multiDens <- function(size_list, plotTit="", legTxt=NULL, legPos="topright", my_ylab="density", my_xlab="") {
  
  dens <- lapply(size_list, function(x) density(na.omit(x)))
  names(dens) <- names(size_list)
  
  lengthDens <- unlist(lapply(size_list, function(x) length(na.omit(x))))
  
  plot(NA, xlim=range(sapply(dens, "[", "x")), ylim=range(sapply(dens, "[", "y")), 
       main=plotTit, xlab=my_xlab, ylab=my_ylab)
  foo <- mapply(lines, dens, col=1:length(dens))
  if(is.null(legTxt)){
    # legTxt <- names(dens)
    legTxt <- paste0(names(dens), " (n=", lengthDens, ")")
  }
  legend(legPos, legend=legTxt, fill=1:length(dens), bty='n')
}



######################################################################################################################################################################################################
###################################################################################################################################################################################################### calculate_MoC_with_domainTypes (function)
######################################################################################################################################################################################################

calculate_MoC_with_domainTypes <- function(set1, set2, chr_len, gaps_as_clusters = FALSE, file_as_input = FALSE) {
#******************************************************************************************** FUNCTION DEFINITION
    getOverlap <- function(a, b){
        # smallest end - biggest start    
#        cat(paste0("a=",a[1],"/", a[2], "\n"))
#                cat(paste0("b=",b[1],"/", b[2], "\n"))
#        overlap = max(0, min(a[2], b[2]) - max(a[1], b[1]) +1)
#        cat(paste0("overlap=", overlap, "\n"))

        return(max(0, min(a[2], b[2]) - max(a[1], b[1]) +1))
    }
    prepareDT <- function(dt, chr_len) {
        stopifnot(is.numeric(chr_len))
        stopifnot(ncol(dt) == 3)
        # do not take the 1st column with the "chr6"    
        dt <- dt[,2:3]
        colnames(dt) <- c("Start", "End")
        # add a column filled with 0
        dt$is_gap <- 0
        # test that nrow dt is bigger than 1 otherwise the 2:... will create NA
        if(nrow(dt) > 1) {
          # create data frame that will hold the gap for the 1st dataset
          # start of the gap = end of the end of the TAD + 1 (do not take the last row)    
          # end of the gap = start of the TAD - 1 (do not take the first row)
          dt_gaps <- data.frame( Start = (dt$End[1:(nrow(dt)-1)] + 1),
                                  End = (dt$Start[2:nrow(dt)] -1))
          stopifnot(is.numeric(dt_gaps$Start[1]))
          stopifnot(is.numeric(dt_gaps$End[1]))    
          # select only the row with end > start
          dt_gaps <- dt_gaps[dt_gaps$Start < dt_gaps$End,]
        } else{
          dt_gaps <- data.frame(Start = numeric(), End=numeric())
        }
        # ad gaps at the beginning until first TAD and at the end until end of chromo size
        # CHANGE MZ: FIRST DOMAIN APPENDED SHOULD START WITH 1 NOT WITH 0
        #pgaps1 = pgaps1.append(pd.DataFrame([[0, p1.iloc[0,0]-1], [p1.iloc[p1.shape[0]-1,1]+1, chr_len]], columns=['Start', 'End']), ignore_index=True)
        # THERE WAS A PROBLEM IF THE LAST TAD WAS UNTIL THE LAST CHROMO IT WHOULD HAD LAST ROW WITH CHR_LEN+1 CHR_LEN
        #dt_gaps = dt_gaps.append(pd.DataFrame([[1, dt.iloc[0,0]-1], [dt.iloc[dt.shape[0]-1,1]+1, chr_len]], columns=['Start', 'End']), ignore_index=True)
        # if needed, add gap before 1st TAD
        if(dt$Start[1] > 1) {
            tmpDT <- data.frame(Start = 1, End = dt$Start[1]-1)
            dt_gaps <- rbind(tmpDT, dt_gaps)
        }
        # if needed, add gap until chromosome end    
        if(dt$End[nrow(dt)] < chr_len) {
            tmpDT <- data.frame(Start = dt$End[nrow(dt)] + 1, End = chr_len)
            dt_gaps <- rbind(dt_gaps, tmpDT)        
        }
        # add a column to indicate there are gaps
        if(nrow(dt_gaps) > 0) {
            dt_gaps$is_gap <- 1
            dt_final <- rbind(dt, dt_gaps)
        } else{
            dt_final <- dt
        }
        dt_final <- dt_final[order(dt_final$Start),]
        return(dt_final)
    }
#********************************************************************************************
    if(file_as_input) {
        if(file.info(set1)$size == 0 & file.info(set2)$size == 0)
            return(NA)
        if(file.info(set1)$size == 0 & file.info(set2)$size > 0)
            return(0)            
        set1DT <- read.delim(set1, header=F, col.names=c("chromo", "start", "end"))
        set2DT <- read.delim(set2, header=F, col.names=c("chromo", "start", "end"))
    } else {
        set1DT <- set1
        set2DT <- set2
        colnames(set1DT) <- c("chromo", "start", "end")
        colnames(set2DT) <- c("chromo", "start", "end")
    }
    stopifnot(is.numeric(set1DT$start[1]))
    stopifnot(is.numeric(set2DT$start[1]))
    stopifnot(is.numeric(set1DT$end[1]))        
    stopifnot(is.numeric(set2DT$end[1]))            
    
    # prepare data for 1st caller
    ptot1 = prepareDT(set1DT, chr_len)

    # prepare data for 2nd caller
    ptot2 = prepareDT(set2DT, chr_len)

    # number of clusters for each TAD caller (# of rows)
    Nclust1 = nrow(ptot1)
    Nclust2 = nrow(ptot2)

#    cat(paste0("ptot1: ", paste0(dim(ptot1), collapse=","), "\n"))
#    cat(paste0("ptot2: ", paste0(dim(ptot2), collapse=","), "\n"))    

    # in MoC definition: if I = J = 1 => MoC = 1
    if(Nclust1==1 & Nclust2==1)
        return(1)
        
    # added for the problem of negative value if compared 1 single whole-chromosome domain to more than 1 domain
    # in this case return NA
    if( (Nclust1 == 1 & Nclust2 > 1) | (Nclust1 > 1 & Nclust2 == 1) )
        return(NA)

    # interTAD regions are considered as domains
    if(gaps_as_clusters){
        crossum = 0
        # iterate over dt1 -> for each domain retrieve the overlapping domains ("fragments") and compute mutual concordance
        for(i in 1:Nclust1) {
            # index of the 1st end in dt2 that has end > than current start
            #firstj = np.argmin(ptot2["End"][ptot2["End"]>ptot1.iloc[i,0] ])        
            #firstj = which.min(ptot2$End[ptot2$End > ptot1$Start[i]])
            firstj = min(which(ptot2$End > ptot1$Start[i]))
            # index of the 1st start in dt2 with start < current end
            # lastj = np.argmax(ptot2["Start"][ptot2["Start"]<ptot1.iloc[i,1] ])
            # lastj = which.max(ptot2$Start[ptot2$Start < ptot1$End[i]])            
            lastj = max(which(ptot2$Start < ptot1$End[i]))
            
            if(is.infinite(firstj) | is.infinite(lastj)) next
            
            for(j in (firstj:lastj)) {
                crossum = crossum + (getOverlap(ptot1[i,1:2], ptot2[j,1:2]))**2/
                                        ( (ptot1[i,2] - ptot1[i,1] + 1)*(ptot2[j,2] - ptot2[j,1] + 1) )
            }            
        }            
        MoC = 1/(sqrt(Nclust1*Nclust2)-1)*(crossum - 1)
        return(MoC)
    } else{          
    # if not interTAD_as_domains => concordance term for <TAD <-> interTAD regions> = 0
        crossum = 0
        for(i in 1:Nclust1 ){
            #firstj = np.argmin(ptot2["End"][ptot2["End"]>ptot1.iloc[i,0] ])                
#            firstj = which.min(ptot2$End[ptot2$End > ptot1$Start[i]])
            firstj = min(which(ptot2$End > ptot1$Start[i]))
            # lastj = np.argmax(ptot2["Start"][ptot2["Start"]<ptot1.iloc[i,1] ])            
#            lastj = which.max(ptot2$Start[ptot2$Start < ptot1$End[i] ])
            lastj = max(which(ptot2$Start < ptot1$End[i] ))            
            #cat(paste0("firstj = ", firstj, "\n"))
            #cat(paste0("lastj = ", lastj, "\n"))            
            for(j in (firstj:lastj)){
#                    cat(paste0("crossum = " ,crossum,"\n"))            
                # ! xor => 0 if this is not the same (ptot1["is_gap"][i] == ptot2["is_gap"][j])
    #            crossum += ( ~np.logical_xor(ptot1$is_gap[i],ptot2$is_gap[j]) )*(getOverlap(ptot1[i,1:2], ptot2[j,1:2]))**2/ \
    #                           ( (ptot1[i,2] - ptot1[i,1] + 1)*(ptot2[j,2] - ptot2[j,1] + 1) )
                crossum = crossum + ( as.numeric(ptot1$is_gap[i] == ptot2$is_gap[j]) )*(getOverlap(ptot1[i,1:2], ptot2[j,1:2]))**2/
                               ( (ptot1[i,2] - ptot1[i,1] + 1)*(ptot2[j,2] - ptot2[j,1] + 1) )              
            }             
        }
        MoC = 1/(sqrt(Nclust1*Nclust2)-1)*(crossum - 1)
        return(MoC)
    }
}   
