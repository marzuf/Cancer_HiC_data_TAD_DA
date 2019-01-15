startTime <- Sys.time()

library(foreach)
library(doMC)

options(scipen=100)


cat("> START: compare_resol_vs_otherStats.R\n")
# Rscript compare_resol_vs_otherStats.R

SSHFS <- FALSE
setDir <- ifelse(SSHFS, "/media/electron", "")

source(file.path(setDir, "/mnt/etemp/marie/Dixon2018_integrative_data", "utils_fct.R"))

if(SSHFS) setwd("/media/electron/mnt/etemp/marie/Dixon2018_integrative_data")

registerDoMC(ifelse(SSHFS, 2, 40))

outFold <- file.path("COMPARE_RESOL_VS_OTHERSTATS")
dir.create(outFold, recursive = TRUE)

# logFile <- file.path(outFold, "resol_vs_otherStats_logFile.txt")
# if(!SSHFS) system(paste0("rm -f ", logFile))
# if(SSHFS) logFile <- ""

plotType <- "svg"

myHeight <- ifelse(plotType=="png", 500, 7)
myWidth <- myHeight

source("utils_fct.R")

cexLab <- 1.2
cexAxis <- 1.2

matchDT <- eval(parse(text = load(file.path(
  "CMP_DATASETS_MATCHING",
  "all_match_dt.Rdata"
))))
head(matchDT)
# ds1             ds2 chromo strictMatchRatio looseMatchRatio bdMatchRatio
# 1 pipConsensus breastConsensus   chr1        0.4285714       0.5164835    0.7678959
# 2 pipConsensus breastConsensus  chr10        0.4202128       0.5265957    0.7754237
matchDS <- unique(c(matchDT$ds1, matchDT$ds2) )

nbrDT <- eval(parse(text = load(file.path(
  "CMP_DATASETS_NBRTADS",
  "all_nbr_dt.Rdata"
))))
head(nbrDT)
# ds1 chromo chromoCover nTADs meanSizeTADs
# 1 pipConsensus   chr1   0.8411425   364     576044.0
# 2 pipConsensus  chr10   0.8960732   188     645744.7
nbrDS <- unique(nbrDT$ds1)

resolDT <- eval(parse(text = load(file.path(
  "CMP_DATASETS_RESOL",
  "all_resol_DT.Rdata"
))))
head(resolDT)
# dataset chromo  countSum rowAbove1000 countSum_log10 datasetLabel
# 1  GSE105194_ENCFF027IEO   chr1 314471319       0.9660       8.497581     astroCL1
# 12 GSE105194_ENCFF027IEO   chr2 256565529       0.9840       8.409198     astroCL1
resolDS <- unique(resolDT$datasetLabel)

# compare to resolDT, chromo_resolDT holds 1 row per Hi-C matrix row
# chromo_resolDT <- eval(parse(text = load(file.path(
#   "CMP_DATASETS_RESOL",
#   "all_chromo_DT.Rdata"
# ))))
# head(chromo_resolDT)
# # dataset chromo matrixDim rowIdx   rowSum rowSum_log10 rowSumNoOut rowSumNoOut_log10 dataset_label
# # 1 GSE105194_ENCFF027IEO   chr1      6080      1 53554.14     4.728793    53554.14          4.728793      astroCL1
# # 2 GSE105194_ENCFF027IEO   chr1      6080      2 53947.06     4.731968    53947.06          4.731968      astroCL1
# chromo_resolDS <- unique(chromo_resolDT$dataset_label)

cat("setdiff(matchDS, resolDS)\n")
setdiff(matchDS, resolDS)
# [1] "pipConsensus"    "breastConsensus" "mcf7Consensus"   "lungConsensus"   "kidneyConsensus" "skinConsensus"   "astroConsensus" 
cat("setdiff(matchDS, nbrDS)\n")
setdiff(matchDS, nbrDS)
# 
cat("setdiff(resolDS, nbrDS)\n")
setdiff(resolDS, nbrDS)
# [1] "astroCL2" "coloCL2"  "astroCL4" "lympho1" 


# stopifnot(length(setdiff(matchDS, resolDS)) == 0)  # !!! will need to be true at some point !!!
# stopifnot(length(setdiff(matchDS, nbrDS)) == 0)  # !!! will need to be true at some point !!!
# stopifnot(length(setdiff(nbrDS, resolDS)) == 0)  # !!! will need to be true at some point !!!

commonDS <- intersect(intersect(resolDS, nbrDS), matchDS)

agg_resolDT <- resolDT
agg_resolDT$dataset <- agg_resolDT$chromo <- NULL
colnames(agg_resolDT)[colnames(agg_resolDT) == "datasetLabel"] <- "dataset"
agg_resolDT <- aggregate(. ~ dataset, data = agg_resolDT, FUN=mean, na.rm=T)
head(agg_resolDT)
# dataset  countSum rowAbove1000 countSum_log10
# 1  astroCL1 138927507    0.8619217       7.947055
# 2  astroCL2  23945826    0.8686489       7.197720
stopifnot(commonDS %in% agg_resolDT$dataset)
agg_resolDT <- agg_resolDT[agg_resolDT$dataset %in% commonDS,]

agg_nbrDT <- nbrDT 
agg_nbrDT$chromo <- NULL
colnames(agg_nbrDT)[colnames(agg_nbrDT) == "ds1"] <- "dataset"
agg_nbrDT <- aggregate(. ~ dataset, data = agg_nbrDT, FUN=mean, na.rm=T)
head(agg_nbrDT)
# dataset chromoCover    nTADs meanSizeTADs
# 1       astroCL1   0.8498755 221.9565     478902.9
# 2       astroCL3   0.8544148 203.1739     527876.1
# 3 astroConsensus   0.6812792 167.6087     510335.6
stopifnot(commonDS %in% agg_nbrDT$dataset)
agg_nbrDT <- agg_nbrDT[agg_nbrDT$dataset %in% commonDS,]

agg_matchDT_tmp <- matchDT
agg_matchDT_tmp$chromo <- NULL
agg_matchDT_1 <- agg_matchDT_tmp
agg_matchDT_2 <- agg_matchDT_tmp
agg_matchDT_1$dataset <- agg_matchDT_tmp$ds1
agg_matchDT_2$dataset <- agg_matchDT_tmp$ds2
agg_matchDT <- rbind(agg_matchDT_1, agg_matchDT_2)
agg_matchDT$ds1 <- agg_matchDT$ds2 <- NULL
agg_matchDT <- aggregate(. ~ dataset, data = agg_matchDT, FUN=mean, na.rm=T)
head(agg_matchDT)
# dataset strictMatchRatio looseMatchRatio bdMatchRatio
# 1       astroCL1        0.4408448       0.5193087    0.6893682
# 2       astroCL3        0.4319568       0.5052715    0.6809180
stopifnot(commonDS %in% agg_matchDT$dataset)
agg_matchDT <- agg_matchDT[agg_matchDT$dataset %in% commonDS,]

all_stats_DT <- merge(agg_nbrDT, 
                      merge(agg_matchDT, agg_resolDT, by="dataset", all.x=TRUE, all.y=TRUE), 
                      by="dataset", all.x=TRUE, all.y=TRUE)
head(all_stats_DT)
# dataset chromoCover    nTADs meanSizeTADs strictMatchRatio looseMatchRatio bdMatchRatio  countSum rowAbove1000 countSum_log10
# 1  astroCL1   0.8498755 221.9565     478902.9        0.4408448       0.5193087    0.6893682 138927507    0.8619217       7.947055
# 2  astroCL3   0.8544148 203.1739     527876.1        0.4319568       0.5052715    0.6809180 119946585    0.8629870       7.886118
stopifnot(!any(is.na(all_stats_DT)))

#  "dataset"          "chromoCover"      "nTADs"            "meanSizeTADs"    
# "strictMatchRatio" "looseMatchRatio"  "bdMatchRatio"    
# "countSum"         "rowAbove1000"     "countSum_log10"  


resol_vars <- c( "countSum",  
                 "countSum_log10" ,
              "rowAbove1000" )

other_vars <- colnames(all_stats_DT)[!colnames(all_stats_DT) %in% c("dataset", resol_vars)]

all_stats_DT$tissues <- as.factor(gsub("(.+)CL.*", "\\1", all_stats_DT$dataset))

var_names <- c(
                   "chromoCover" = "chromosome coverage",
                   "nTADs" = "# of TADs",
                   "meanSizeTADs" = "mean TAD size",
                    "strictMatchRatio"  = "strict matching ratio",
                    "looseMatchRatio" = "loose matching ratio",
                    "bdMatchRatio" = "boundary matching ratio",
                    "countSum" = "sum of contact counts",
                    "rowAbove1000" = "ratio of rows with > 1000 contact counts",
                    "countSum_log10" = "sum of contact counts [log10]"  )



nDS <- length(commonDS)

resolVar = resol_vars[1]
otherVar = other_vars[1]
all_corr_dt <- foreach(resolVar = resol_vars, .combine='rbind') %do% {
  resolvar_dt <- foreach(otherVar = other_vars, .combine='rbind') %dopar% {
    
    leg_pos <- ifelse(otherVar %in% c("chromoCover", "meanSizeTADs", "nTADs"), "bottomleft", "bottomright")
    
    myx <- all_stats_DT[,resolVar]
    myy <- all_stats_DT[,otherVar]
    
    stopifnot(resolVar %in% names(var_names))
    stopifnot(otherVar %in% names(var_names))
    
    myxlab <- paste0(var_names[resolVar])
    myylab <- paste0(var_names[otherVar])
    
    myTit <- paste0(myylab, " vs. ", myxlab)
    
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
             main = myTit,
            col = as.numeric(all_stats_DT[,"tissues"]),
       cex.axis = cexAxis, cex.lab = cexLab
    )
    text(x=myx,
         y=myy,
         labels=all_stats_DT[,"dataset"], 
         col = as.numeric(all_stats_DT[,"tissues"]),
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
    
    corr_coeff <- as.numeric(cor.test(myx,myy)$estimate)
    stopifnot(!is.na(corr_coeff))
    corr_pval <- as.numeric(cor.test(myx, myy)$p.val)
    stopifnot(!is.na(corr_pval))
    
    data.frame(
      resol_var = resolVar,
      other_var = otherVar,
      PCC = corr_coeff,
      pval = corr_pval,
      stringsAsFactors=FALSE
      )
  } # end iterating over otherVar
  resolvar_dt
} #  end iterating over resolVar


outFile <- file.path(outFold, "all_corr_dt.Rdata")
save(all_corr_dt, file = outFile)
cat(paste0("... written: ", outFile, "\n"))

######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))








