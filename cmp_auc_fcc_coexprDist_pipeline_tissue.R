startTime <- Sys.time()

library(foreach)
library(doMC)
library(ggplot2)
library(stringi)

options(scipen=100)

cat("> START: cmp_auc_fcc_coexprDist_pipeline_tissue.R\n")
# Rscript cmp_auc_fcc_coexprDist_pipeline_tissue.R

SSHFS <- FALSE
setDir <- ifelse(SSHFS, "/media/electron", "")

coexprdistFolder <- file.path("AUC_COEXPRDIST_WITHFAM_SORTNODUP")
stopifnot(dir.exists(coexprdistFolder))
# AUC_COEXPRDIST_WITHFAM_SORTNODUP/GSM1631185_MCF7_vs_GSE75070_MCF7_shGFP/TCGAbrca_lum_bas_hgnc/hgnc_family_short/auc_values.Rdata

coexprdistOldFolder <-  file.path(
  setDir,
  "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom",
  "AUC_COEXPRDIST_WITHFAM_SORTNODUP"
)
stopifnot(dir.exists(coexprdistOldFolder))


old_fcc_folder <- file.path(
  setDir,
  "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom",
  "OUTPUT_FOLDER"
)
stopifnot(dir.exists(old_fcc_folder))

fcc_folder <- file.path("PIPELINE", "OUTPUT_FOLDER")
stopifnot(dir.exists(fcc_folder))

plotType <- "png"
myHeight <- ifelse(plotType == "png", 500, 7)
myWidth <- myHeight
myGGheight <- ifelse(plotType == "png", 500, 7)
myGGwidth <- ifelse(plotType == "png", 700, 10)
plotCex <- 1.2

outFold <- file.path("CMP_AUC_FCC_COEXPRDIST_PIPELINE_TISSUE")
dir.create(outFold, recursive = TRUE)

source( file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg"), "set_dataset_colors.R"))
head(score_DT)

source( file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_11_18"), "analysis_utils.R"))

dataset_proc_colors <- setNames(score_DT$proc_col, score_DT$dataset)
length(dataset_proc_colors)

my_colors_leg <- my_colors

#### BUILD THE BARPLOT AUC FCC RANKING 

# retrieve all the FCC auc ratios available

all_ratio_files <- list.files(coexprdistFolder,
                          pattern = "auc_values.Rdata",
                          recursive=TRUE,
                          full.names=TRUE)
stopifnot(length(all_ratio_files) > 0)

curr_file = "AUC_COEXPRDIST_WITHFAM_SORTNODUP/GSM1631185_MCF7_vs_GSE75070_MCF7_shGFP/TCGAbrca_lum_bas_hgnc/hgnc_family_short/auc_values.Rdata"

all_auc_DT <- foreach(curr_file = all_ratio_files, .combine='rbind') %dopar% {
  curr_ds <- basename(dirname(dirname(curr_file)))
  ds_name <- gsub("(.+?_.+?_.+?)_.+", "\\1", curr_ds)
  stopifnot(grepl("TCGA", curr_ds))
  all_ratios <- eval(parse(text = load(curr_file)))
  aucCoexprDist <- as.numeric(all_ratios["auc_ratio_same_over_diff_distVect"])
  stopifnot(!is.na(aucCoexprDist))
  aucCoexprDistSameFam <- as.numeric(all_ratios["auc_ratio_sameFam_same_over_diff_distVect"])
  stopifnot(!is.na(aucCoexprDistSameFam))
  
  old_file <- list.files(file.path(coexprdistOldFolder,
                        curr_ds), recursive = TRUE, full.names = TRUE,
                        pattern = "auc_values.Rdata")
  stopifnot(length(old_file) == 1)
  stopifnot(file.exists(old_file))
  
  old_ratios <- eval(parse(text = load(old_file)))
  old_aucCoexprDist <- as.numeric(old_ratios["auc_ratio_same_over_diff_distVect"])
  old_aucCoexprDistSameFam <- as.numeric(old_ratios["auc_ratio_sameFam_same_over_diff_distVect"])
  stopifnot(! is.na(old_aucCoexprDist))
  stopifnot(! is.na(old_aucCoexprDistSameFam))
  
  
  curr_fcc_file <- list.files(file.path(fcc_folder),
                                pattern = paste0("auc_ratios.Rdata"),
                                recursive=TRUE,
                                full.names=TRUE)
  curr_fcc_file <- curr_fcc_file[grepl(ds_name, curr_fcc_file)]
  stopifnot(length(curr_fcc_file) == 1)
  stopifnot(file.exists(curr_fcc_file))
  
  all_fcc_ratios <- eval(parse(text = load(curr_fcc_file)))
  aucfcc <- as.numeric(all_fcc_ratios["prodSignedRatio_auc_permGenes"])
  stopifnot(! is.na(aucfcc))
  
  
  old_fcc_file <- file.path(old_fcc_folder,
                        ds_name,
                        "170revision2EZH2_score_auc_pval_permGenes",
                        "auc_ratios.Rdata")
  stopifnot(file.exists(old_fcc_file))
  old_ratios <- eval(parse(text = load(old_fcc_file)))
  old_aucfcc <- as.numeric(old_ratios["prodSignedRatio_auc_permGenes"])
  stopifnot(! is.na(old_aucfcc))
  
  
  data.frame(
    dataset = ds_name,
    aucFCC = aucfcc,
    old_aucFCC = old_aucfcc,
    aucCoexprDist = aucCoexprDist,
    old_aucCoexprDist = old_aucCoexprDist,
    aucCoexprDistSameFam = aucCoexprDistSameFam,
    old_aucCoexprDistSameFam = old_aucCoexprDistSameFam, 
    stringsAsFactors = FALSE
  )
}


#################################### 
### correlation FCC ~ COEXPR DIST -- OLD SCORES
#################################### 

all_comps <- list(

  c("aucFCC", "old_aucFCC"),
  c("aucCoexprDist", "old_aucCoexprDist"),
  c("aucCoexprDistSameFam", "old_aucCoexprDistSameFam"),

  c("aucFCC", "aucCoexprDist"),
  c("aucFCC", "aucCoexprDistSameFam"),
  c("old_aucFCC", "old_aucCoexprDist"),
  c("old_aucFCC", "old_aucCoexprDistSameFam")

)

i=1
for( i in seq_along(all_comps) ) {
  
  xvar <- all_comps[[i]][1]
  yvar <- all_comps[[i]][2]
  
  myx <- 100*(all_auc_DT[,xvar]-1)
  myy <- 100*(all_auc_DT[,yvar]-1)
  
  leg_pos <- "topleft"
  
  if(grepl("SameFam", xvar) & grepl("SameFam", yvar))
    leg_pos <- "topright"
  
  if(grepl("CoexprDist$", xvar) & grepl("CoexprDist$", yvar))
    leg_pos <- "topright"
  
  if(xvar=="aucFCC" & yvar=="aucCoexprDist")
    leg_pos <- "topright"
  
#  if(grepl("old_", xvar)) {
#    stopifnot(grepl("old_", yvar))
#    labpart <- "pipeline TADs"
#    xvar <- gsub("old_", "", xvar)
#    yvar <- gsub("old_", "", yvar)
#    
#  } else {
#    labpart <- "tissue TADs"
#  }
  if(grepl("old_", xvar)) {
    labpartX <- "pipeline TADs"
    xvar <- gsub("old_", "", xvar)
  } else {
    labpartX <- "tissue TADs"
  }
  if(grepl("old_", yvar)) {
    labpartY <- "pipeline TADs"
    yvar <- gsub("old_", "", yvar)
  } else {
    labpartY <- "tissue TADs"
  }
 
  
#  myTit <- paste0(" % increase AUC\n", xvar, " vs. ", yvar, " (", gsub(" ", "", labpart), ")")
  myTit <- paste0(" % increase AUC\n", xvar, "  (", gsub(" ", "", labpartX), ") vs. ", yvar, " (", gsub(" ", "", labpartY), ")")
  #myTit <- paste0("Tissue-specific consensus vs. pipeline consensus % increase FCC AUC")
#  x_lab <- paste("% increase", xvar, " (", labpart, ")")
#  y_lab <- paste("% increase", yvar, " (", labpart, ")")
  x_lab <- paste("% increase", xvar, " (", labpartX, ")")
  y_lab <- paste("% increase", yvar, " (", labpartY, ")")
  mySub <- paste0("(# DS = ", nrow(all_auc_DT), ")")
  
  mynames <- all_auc_DT$dataset
  
  curr_colors <- as.character(cancer_subColors[as.character(cancer_subAnnot[mynames])])
  stopifnot(!is.na(curr_colors))
  
#  outFile <- file.path(outFold, paste0(xvar, "_", yvar, "_", gsub(" ", "", labpart), ".", plotType))
  outFile <- file.path(outFold, paste0(xvar, "_", gsub(" ", "", labpartX),"_", yvar, "_", gsub(" ", "", labpartY), ".", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  plot(x = myx,
       y = myy,
       main = myTit,
       xlab=paste0(x_lab),
       ylab=paste0(y_lab),
       cex = 0.7, pch=16,
       cex.lab = plotCex, cex.axis = plotCex
  )
  mtext(side=3, text = mySub)
  text(x = myx,
       y = myy,
       labels = mynames,
       col = curr_colors,
       pos=3, cex = 0.7)
  addCorr(x=myx, 
          y=myy,
          legPos="bottomright", 
          corMet="spearman",
          bty="n") 
  legend(leg_pos,
         legend=unique(cancer_subAnnot[mynames]), #names(curr_colors),
         lty=1,
         col = unique(curr_colors),
         lwd = 5,
         bty="n",
         cex = 0.7)
  add_curv_fit(x = myx, 
               y=myy,
               withR2 = FALSE, R2shiftX = -0.03, R2shiftY = 0, col="grey", lty=2)
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
}


  
