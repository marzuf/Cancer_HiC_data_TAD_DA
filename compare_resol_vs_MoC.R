startTime <- Sys.time()

library(foreach)
library(doMC)
library(dplyr)

options(scipen=100)

cat("> START: compare_resol_vs_MoC.R\n")
# Rscript compare_resol_vs_MoC.R

SSHFS <- FALSE
setDir <- ifelse(SSHFS, "/media/electron", "")
setDir <- ifelse(SSHFS, "~/media/electron", "")

source(file.path("utils_fct.R"))

# if(SSHFS) setwd("/media/electron/mnt/etemp/marie/Cancer_HiC_data_TAD_DA")
if(SSHFS) setwd("~/media/electron/mnt/etemp/marie/Cancer_HiC_data_TAD_DA")

registerDoMC(ifelse(SSHFS, 2, 40))

mocFile <- file.path(
  "CMP_DATASETS_MOC",
  "all_MoC_dt.Rdata"
)
stopifnot(file.exists(mocFile))
resolFile <- file.path(
  "CMP_DATASETS_RESOL",
  "all_resol_DT.Rdata"
)
stopifnot(file.exists(resolFile))

outFold <- file.path("COMPARE_RESOL_VS_MOC")
dir.create(outFold, recursive = TRUE)

# logFile <- file.path(outFold, "resol_vs_otherStats_logFile.txt")
# if(!SSHFS) system(paste0("rm -f ", logFile))
# if(SSHFS) logFile <- ""

plotType <- "svg"

myHeight <- ifelse(plotType=="png", 500, 7)
myWidth <- myHeight

mocDT <- eval(parse(text = load(mocFile)))
head(mocDT)
# ds1             ds2 chromo       MoC
# 1 pipConsensus breastConsensus   chr1 0.6174482
# 2 pipConsensus breastConsensus  chr10 0.6232023
# 3 pipConsensus breastConsensus  chr11 0.5779753
mocDS <- unique(c(mocDT$ds1, mocDT$ds2) )

resolDT <- eval(parse(text = load(resolFile)))
head(resolDT)
# dataset chromo  countSum rowAbove1000 countSum_log10 datasetLabel
# 1  GSE105194_ENCFF027IEO   chr1 314471319       0.9660       8.497581     astroCL1
# 12 GSE105194_ENCFF027IEO   chr2 256565529       0.9840       8.409198     astroCL1
resolDS <- unique(resolDT$datasetLabel)

cat("setdiff(mocDS, resolDS)\n")
setdiff(mocDS, resolDS)

commonDS <- intersect(mocDS, resolDS)

mocDT <- mocDT[mocDT$ds1 %in% commonDS & mocDT$ds2 %in% commonDS,]
resolDT <- resolDT[resolDT$datasetLabel %in% commonDS,]

### HAVE i THE MOC FOR EACH CHROMO IN MOC DT ???

resolDT$chromo <- as.character(resolDT$chromo)
mocDT$chromo <- as.character(mocDT$chromo)

resol1DT <- resolDT
resol1DT$dataset <- NULL
colnames(resol1DT) <- paste0(colnames(resol1DT), "_ds1")
merge1 <- left_join(mocDT, resol1DT, 
                    by = c("ds1" = "datasetLabel_ds1", "chromo" ="chromo_ds1"))  # or datasetLabel_ds1 ??? if chromo available in MoC dt, otherwise just dataset

resol2DT <- resolDT
resol2DT$dataset <- NULL
colnames(resol2DT) <- paste0(colnames(resol2DT), "_ds2")
merge2 <- left_join(merge1, resol2DT, by = c("ds2" = "datasetLabel_ds2", "chromo" ="chromo_ds2"))  # or datasetLabel_ds1 ??? if chromo available in MoC dt, otherwise just dataset

resol_MoC_DT <- merge2

mycolumns <- c("rowAbove1000", "countSum", "countSum_log10")

otherVar <- "MoC"


var_names <- c(
                    "MoC" = "MoC",
                    "countSum_ratio" = "ratio of sum of contact counts",
                    "rowAbove1000_ratio" = "ratio of ratio of rows with > 1000 contact counts",
                    "countSum_log10_ratio" = "ratio sum of contact counts [log10]"  
)

# here:2datasets,not1signletissue
# resol_MoC_DT$tissues <- as.factor(gsub("(.+)CL.*", "\\1", resol_MoC_DT$ds1))

nDS <- length(commonDS)

resol_MoC_DT$dataset <- paste0(resol_MoC_DT$ds1, " vs. ", resol_MoC_DT$ds2)

### 1 row for each chromo
for(mycol in mycolumns) {
     resolVar <- paste0(mycol, "_ratio")
     resol_MoC_DT[,resolVar] <- resol_MoC_DT[,paste0(mycol, "_ds1")] / resol_MoC_DT[,paste0(mycol, "_ds2")]
    # leg_pos <- ifelse(mycol %in% c("TOSET"), "bottomleft", "bottomright")
    leg_pos <- "topright"
    myx <- resol_MoC_DT[,resolVar]
    myy <- resol_MoC_DT[,otherVar]
    stopifnot(resolVar %in% names(var_names))
    stopifnot(otherVar %in% names(var_names))
    myxlab <- paste0(var_names[resolVar])
    myylab <- paste0(var_names[otherVar])
    myTit <- paste0(myylab, " vs. ", myxlab," (all chromo)")
    mySub <- paste0("(# DS = ", nDS, ")")
    outFile <- file.path(outFold, paste0(otherVar, "_vs_", resolVar, ".", plotType))
    do.call(plotType, list(outFile, height=myHeight, width=myWidth))
    # densplot(
      plot(
            x=myx,
             y=myy,
             xlab=myxlab,
             ylab=myylab,
             pch = 16, cex = 0.7,
            # col = as.numeric(resol_MoC_DT[,"tissues"]),
             main = myTit
    )
    text(x=myx,
         y=myy,
         labels=resol_MoC_DT[,"dataset"], 
         # col = as.numeric(resol_MoC_DT[,"tissues"]),
         cex=0.7)  
    mtext(side=3, text = mySub)
    add_curv_fit(x = myx,
                 y = myy,
                 withR2 = FALSE, lty=2, col="darkgray")
    addCorr(x=myx,
            y=myy,
            bty="n",
            legPos=leg_pos)
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
}

### aggregate -> 1 row for each dataset pair
meanDS_resol_MoC_DT <- resol_MoC_DT
meanDS_resol_MoC_DT$ds1 <- meanDS_resol_MoC_DT$ds2 <- meanDS_resol_MoC_DT$chromo <- NULL
meanDS_resol_MoC_DT <- aggregate(. ~ dataset, data = meanDS_resol_MoC_DT, FUN=mean, na.rm=T)

### 1 row for each chromo
for(mycol in mycolumns) {
  resolVar <- paste0(mycol, "_ratio")
  meanDS_resol_MoC_DT[,resolVar] <- meanDS_resol_MoC_DT[,paste0(mycol, "_ds1")] / meanDS_resol_MoC_DT[,paste0(mycol, "_ds2")]
  # leg_pos <- ifelse(mycol %in% c("TOSET"), "bottomleft", "bottomright")
  leg_pos <- "bottomright"
  myx <- meanDS_resol_MoC_DT[,resolVar]
  myy <- meanDS_resol_MoC_DT[,otherVar]
  stopifnot(resolVar %in% names(var_names))
  stopifnot(otherVar %in% names(var_names))
  myxlab <- paste0(var_names[resolVar])
  myylab <- paste0(var_names[otherVar])
  myTit <- paste0(myylab, " vs. ", myxlab, " (avg. chromo)")
  mySub <- paste0("(# DS = ", nDS, ")")
  outFile <- file.path(outFold, paste0(otherVar, "_vs_", resolVar, "_avgChromo.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  # densplot(
  plot(
    x=myx,
    y=myy,
    xlab=myxlab,
    ylab=myylab,
    pch = 16, cex = 0.7,
    # col = as.numeric(meanDS_resol_MoC_DT[,"tissues"]),
    main = myTit
  )
  text(x=myx,
       y=myy,
       labels=meanDS_resol_MoC_DT[,"dataset"], 
       # col = as.numeric(meanDS_resol_MoC_DT[,"tissues"]),
       cex=0.7)  
  mtext(side=3, text = mySub)
  add_curv_fit(x = myx,
               y = myy,
               withR2 = FALSE, lty=2, col="darkgray")
  addCorr(x=myx,
          y=myy,
          bty="n",
          legPos=leg_pos)
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
}


## stopifnot(length(setdiff(mocDS, resolDS)) == 0)  # !!! will need to be true at some point !!!

#commonDS <- intersect(resolDS, mocDS)

#agg_resolDT <- resolDT
#agg_resolDT$dataset <- agg_resolDT$chromo <- NULL
#colnames(agg_resolDT)[colnames(agg_resolDT) == "datasetLabel"] <- "dataset"
#agg_resolDT <- aggregate(. ~ dataset, data = agg_resolDT, FUN=mean, na.rm=T)
#head(agg_resolDT)
## dataset  countSum rowAbove1000 countSum_log10
## 1  astroCL1 138927507    0.8619217       7.947055
## 2  astroCL2  23945826    0.8686489       7.197720
#stopifnot(commonDS %in% agg_resolDT$dataset)
#agg_resolDT <- agg_resolDT[agg_resolDT$dataset %in% commonDS,]

#agg_nbrDT <- nbrDT 
#agg_nbrDT$chromo <- NULL
#colnames(agg_nbrDT)[colnames(agg_nbrDT) == "ds1"] <- "dataset"
#agg_nbrDT <- aggregate(. ~ dataset, data = agg_nbrDT, FUN=mean, na.rm=T)
#head(agg_nbrDT)
## dataset chromoCover    nTADs meanSizeTADs
## 1       astroCL1   0.8498755 221.9565     478902.9
## 2       astroCL3   0.8544148 203.1739     527876.1
## 3 astroConsensus   0.6812792 167.6087     510335.6
#stopifnot(commonDS %in% agg_nbrDT$dataset)
#agg_nbrDT <- agg_nbrDT[agg_nbrDT$dataset %in% commonDS,]

#agg_mocDT_tmp <- mocDT
#agg_mocDT_tmp$chromo <- NULL
#agg_mocDT_1 <- agg_mocDT_tmp
#agg_mocDT_2 <- agg_mocDT_tmp
#agg_mocDT_1$dataset <- agg_mocDT_tmp$ds1
#agg_mocDT_2$dataset <- agg_mocDT_tmp$ds2
#agg_mocDT <- rbind(agg_mocDT_1, agg_mocDT_2)
#agg_mocDT$ds1 <- agg_mocDT$ds2 <- NULL
#agg_mocDT <- aggregate(. ~ dataset, data = agg_mocDT, FUN=mean, na.rm=T)
#head(agg_mocDT)
## dataset strictMatchRatio looseMatchRatio bdMatchRatio
## 1       astroCL1        0.4408448       0.5193087    0.6893682
## 2       astroCL3        0.4319568       0.5052715    0.6809180
#stopifnot(commonDS %in% agg_mocDT$dataset)
#agg_mocDT <- agg_mocDT[agg_mocDT$dataset %in% commonDS,]

#resol_MoC_DT <- merge(agg_nbrDT, 
#                      merge(agg_mocDT, agg_resolDT, by="dataset", all.x=TRUE, all.y=TRUE), 
#                      by="dataset", all.x=TRUE, all.y=TRUE)
#head(resol_MoC_DT)
## dataset chromoCover    nTADs meanSizeTADs strictMatchRatio looseMatchRatio bdMatchRatio  countSum rowAbove1000 countSum_log10
## 1  astroCL1   0.8498755 221.9565     478902.9        0.4408448       0.5193087    0.6893682 138927507    0.8619217       7.947055
## 2  astroCL3   0.8544148 203.1739     527876.1        0.4319568       0.5052715    0.6809180 119946585    0.8629870       7.886118
#stopifnot(!any(is.na(resol_MoC_DT)))

##  "dataset"          "chromoCover"      "nTADs"            "meanSizeTADs"    
## "strictMatchRatio" "looseMatchRatio"  "bdMatchRatio"    
## "countSum"         "rowAbove1000"     "countSum_log10"  


#resol_vars <- c( "countSum",  
#                 "countSum_log10" ,
#              "rowAbove1000" )

#other_vars <- colnames(resol_MoC_DT)[!colnames(resol_MoC_DT) %in% c("dataset", resol_vars)]

#resol_MoC_DT$tissues <- as.factor(gsub("(.+)CL.*", "\\1", resol_MoC_DT$dataset))

#var_names <- c(
#                   "chromoCover" = "chromosome coverage",
#                   "nTADs" = "# of TADs",
#                   "meanSizeTADs" = "mean TAD size",
#                    "strictMatchRatio"  = "strict matching ratio",
#                    "looseMatchRatio" = "loose matching ratio",
#                    "bdMatchRatio" = "boundary matching ratio",
#                    "countSum" = "sum of contact counts",
#                    "rowAbove1000" = "ratio of rows with > 1000 contact counts",
#                    "countSum_log10" = "sum of contact counts [log10]"  )


#resol_MoC_DT$tissues <- as.factor(gsub("(.+)CL.*", "\\1", resol_MoC_DT$dataset))
#nDS <- length(commonDS)

#resolVar = resol_vars[1]
#otherVar = other_vars[1]
#all_corr_dt <- foreach(resolVar = resol_vars, .combine='rbind') %do% {
#  resolvar_dt <- foreach(otherVar = other_vars, .combine='rbind') %dopar% {
#    
#    leg_pos <- ifelse(otherVar %in% c("chromoCover", "meanSizeTADs", "nTADs"), "bottomleft", "bottomright")
#    
#    myx <- resol_MoC_DT[,resolVar]
#    myy <- resol_MoC_DT[,otherVar]
#    
#    stopifnot(resolVar %in% names(var_names))
#    stopifnot(otherVar %in% names(var_names))
#    
#    myxlab <- paste0(var_names[resolVar])
#    myylab <- paste0(var_names[otherVar])
#    
#    myTit <- paste0(myylab, " vs. ", myxlab)
#    
#    mySub <- paste0("(# DS = ", nDS, ")")
#    
#    outFile <- file.path(outFold, paste0(otherVar, "_vs_", resolVar, ".", plotType))
#    do.call(plotType, list(outFile, height=myHeight, width=myWidth))
#    # densplot(
#      plot(
#            x=myx,
#             y=myy,
#             xlab=myxlab,
#             ylab=myylab,
#             pch = 16, cex = 0.7,
#             main = myTit,
#            col = as.numeric(resol_MoC_DT[,"tissues"])
#    )
#    text(x=myx,
#         y=myy,
#         labels=resol_MoC_DT[,"dataset"], 
#         col = as.numeric(resol_MoC_DT[,"tissues"]),
#         cex=0.7)  
#    mtext(side=3, text = mySub)
#    add_curv_fit(x = myx,
#                 y = myy,
#                 withR2 = FALSE, lty=2, col="darkgray")
#    addCorr(x=myx,
#            y=myy,
#            bty="n",
#            legPos=leg_pos)
#    
#    foo <- dev.off()
#    cat(paste0("... written: ", outFile, "\n"))
#    
#    corr_coeff <- as.numeric(cor.test(myx,myy)$estimate)
#    stopifnot(!is.na(corr_coeff))
#    corr_pval <- as.numeric(cor.test(myx, myy)$p.val)
#    stopifnot(!is.na(corr_pval))
#    
#    data.frame(
#      resol_var = resolVar,
#      other_var = otherVar,
#      PCC = corr_coeff,
#      pval = corr_pval,
#      stringsAsFactors=FALSE
#      )
#  } # end iterating over otherVar
#  resolvar_dt
#} #  end iterating over resolVar


#outFile <- file.path(outFold, "all_corr_dt.Rdata")
#save(all_corr_dt, file = outFile)
#cat(paste0("... written: ", outFile, "\n"))

######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))








