startTime <- Sys.time()

library(foreach)
library(doMC)

options(scipen=100)

cat("> START: compare_resol_vs_otherStats.R\n")
# Rscript compare_resol_vs_otherStats.R

SSHFS <- FALSE
setDir <- ifelse(SSHFS, "/media/electron", "")

source(file.path("utils_fct.R"))
source(file.path("datasets_settings.R"))

# if(SSHFS) setwd("/media/electron/mnt/etemp/marie/Cancer_HiC_data_TAD_DA")
if(SSHFS) setwd("~/media/electron/mnt/etemp/marie/Cancer_HiC_data_TAD_DA")

registerDoMC(ifelse(SSHFS, 2, 40))

outFold <- file.path("COMPARE_RESOL_VS_OTHERSTATS")
dir.create(outFold, recursive = TRUE)

# logFile <- file.path(outFold, "resol_vs_otherStats_logFile.txt")
# if(!SSHFS) system(paste0("rm -f ", logFile))
# if(SSHFS) logFile <- ""

plotType <- "svg"
myHeight <- ifelse(plotType=="png", 500, 7)
myWidth <- myHeight

plotCex <- 1.2

matchDT <- eval(parse(text = load(file.path(
  "CMP_DATASETS_MATCHING",
  "all_match_dt.Rdata"
))))
head(matchDT)
# ds1              ds2 chromo strictMatchRatio looseMatchRatio bdMatchRatio
# 1 MCF-7 ENCSR549MGQ_T47D   chr1        0.4365079       0.6216931    0.7985782
# 2 MCF-7 ENCSR549MGQ_T47D  chr10        0.3006536       0.5490196    0.7455621
matchDS <- unique(c(matchDT$ds1, matchDT$ds2) )

nbrDT <- eval(parse(text = load(file.path(
  "CMP_DATASETS_NBRTADS",
  "all_nbr_dt.Rdata"
))))
head(nbrDT)
# ds1 chromo chromoCover nTADs meanSizeTADs
# 1 MCF-7   chr1   0.8913497   378     587724.9
# 2 MCF-7  chr10   0.9598583   153     850196.1
nbrDS <- unique(nbrDT$ds1)

resolDT <- eval(parse(text = load(file.path(
  "CMP_DATASETS_RESOL",
  "all_resol_DT.Rdata"
))))
head(resolDT)
# dataset chromo countSum rowAbove1000 countSum_log10 datasetLabel
# 1  ENCSR079VIJ_G401   chr1 19193942       0.8917       7.283164    kidneyCl1
# 12 ENCSR079VIJ_G401   chr2 22877966       0.9712       7.359417    kidneyCl1
resolDS <- unique(resolDT$dataset)

cat("setdiff(matchDS, resolDS)\n")
setdiff(matchDS, resolDS)
cat("setdiff(matchDS, nbrDS)\n")
setdiff(matchDS, nbrDS)
cat("setdiff(resolDS, nbrDS)\n")
setdiff(resolDS, nbrDS)

commonDS <- intersect(intersect(resolDS, nbrDS), matchDS)
stopifnot(length(commonDS) > 0)

agg_resolDT <- resolDT
agg_resolDT$datasetLabel <- agg_resolDT$chromo <- NULL
colnames(agg_resolDT)[colnames(agg_resolDT) == "dataset"] <- "dataset"
agg_resolDT <- aggregate(. ~ dataset, data = agg_resolDT, FUN=mean, na.rm=T)
head(agg_resolDT)
stopifnot(commonDS %in% agg_resolDT$dataset)
agg_resolDT <- agg_resolDT[agg_resolDT$dataset %in% commonDS,]
stopifnot(nrow(agg_resolDT) > 0)

agg_nbrDT <- nbrDT 
agg_nbrDT$chromo <- NULL
colnames(agg_nbrDT)[colnames(agg_nbrDT) == "ds1"] <- "dataset"
agg_nbrDT <- aggregate(. ~ dataset, data = agg_nbrDT, FUN=mean, na.rm=T)
head(agg_nbrDT)
stopifnot(commonDS %in% agg_nbrDT$dataset)
agg_nbrDT <- agg_nbrDT[agg_nbrDT$dataset %in% commonDS,]
stopifnot(nrow(agg_nbrDT) > 0)

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
stopifnot(commonDS %in% agg_matchDT$dataset)
agg_matchDT <- agg_matchDT[agg_matchDT$dataset %in% commonDS,]
stopifnot(nrow(agg_matchDT) > 0)

all_stats_DT <- merge(agg_nbrDT, 
                      merge(agg_matchDT, agg_resolDT, by="dataset", all.x=TRUE, all.y=TRUE), 
                      by="dataset", all.x=TRUE, all.y=TRUE)
stopifnot(nrow(all_stats_DT) > 0)
head(all_stats_DT)
# dataset chromoCover    nTADs meanSizeTADs strictMatchRatio looseMatchRatio bdMatchRatio countSum rowAbove1000 countSum_log10
# 1     ENCSR079VIJ_G401   0.8922913 228.2174     526880.7        0.4872456       0.5876920    0.7698598 10335754   0.90403913       6.920983
# 2  ENCSR105KFX_SK-N-DZ   0.8946099 196.8261     607169.0        0.3387719       0.4516679    0.6818796  2497151   0.08507826       6.283619
stopifnot(!any(is.na(all_stats_DT)))
# COLUMNS:
#  "dataset"          "chromoCover"      "nTADs"            "meanSizeTADs"    
# "strictMatchRatio" "looseMatchRatio"  "bdMatchRatio"    
# "countSum"         "rowAbove1000"     "countSum_log10"  

all_stats_DT$datasetID <- all_stats_DT$dataset
stopifnot( all_stats_DT$dataset %in% cl_names)
# all_stats_DT$dataset <- sapply(all_stats_DT$datasetID, function(x) {
#     dslab <- as.character(cl_labs[x])
#     stopifnot(length(dslab) == 1)
#     if(is.na(dslab)) {
#       paste0(names(cl_names[cl_names == x]))
#     } else {
#       paste0(names(cl_names[cl_names == x]), "\n(", dslab, ")")
#     }
# })  
all_stats_DT$dataset <- sapply(all_stats_DT$datasetID, function(x) as.character(names(cl_names[cl_names == x])))
  
  
resol_vars <- c( "countSum",  
                 "countSum_log10" ,
              "rowAbove1000" )

other_vars <- colnames(all_stats_DT)[!colnames(all_stats_DT) %in% c("dataset", "datasetID", resol_vars)]

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
       cex.axis = plotCex, cex.lab = plotCex
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








