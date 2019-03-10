startTime <- Sys.time()

library(foreach)
library(doMC)
library(ggplot2)

options(scipen=100)

cat("> START: cmp_auc_coexprDist_tissue.R\n")
# Rscript cmp_auc_coexprDist_tissue.R

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


plotType <- "svg"
myHeight <- ifelse(plotType == "png", 300, 7)
myWidth <- myHeight
myGGheight <- ifelse(plotType == "png", 300, 7)
myGGwidth <- ifelse(plotType == "png", 500, 10)
plotCex <- 1.2

outFold <- file.path("CMP_AUC_COEXPRDIST_TISSUE")
dir.create(outFold, recursive = TRUE)

source( file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg"), "set_dataset_colors.R"))
head(score_DT)

source("utils_fct.R")
source("colors_utils.R")

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

curr_file = "AUC_COEXPRDIST_WITHFAM_SORTNODUP/ENCSR079VIJ_G401_40kb/TCGAkich_norm_kich_hgnc/hgnc_family_short/auc_values.Rdata"

all_aucCoexprDist_DT <- foreach(curr_file = all_ratio_files, .combine='rbind') %dopar% {
  curr_ds <- basename(dirname(dirname(curr_file)))
  ds_name <- gsub("(.+?_.+?_.+?)_.+", "\\1", curr_ds)
  stopifnot(grepl("TCGA", curr_ds))
  all_ratios <- eval(parse(text = load(curr_file)))
  aucCoexprDist <- as.numeric(all_ratios["auc_ratio_same_over_diff_distVect"])
  stopifnot(!is.na(aucCoexprDist))
  aucCoexprDistSameFam <- as.numeric(all_ratios["auc_ratio_sameFam_same_over_diff_distVect"])
  stopifnot(!is.na(aucCoexprDistSameFam))
  


    
  data.frame(
    dataset = ds_name,
    datasetLabel = paste0(basename(dirname(dirname(dirname(curr_file)))), "\n", basename(dirname(dirname(curr_file)))),
    aucCoexprDist = aucCoexprDist,
    aucCoexprDistSameFam = aucCoexprDistSameFam,
    stringsAsFactors = FALSE
  )
}








##############################################################################################
############################################################################################## BARPLOT WITH ONE AUC RATIOS
##############################################################################################

var_names <- c(aucCoexprDist = "% increase AUC coexprDist (tissue-specific)",
               old_aucCoexprDist = "% increase AUC coexprDist (pipeline consensus)",
               aucCoexprDistSameFam = "% increase AUC coexprDist (same fam.; tissue-specific)",
               old_aucCoexprDistSameFam = "% increase AUC coexprDist (same fam.; pipeline consensus)"
               )

coexprdist="old_aucCoexprDist"

for(coexprdist in c("aucCoexprDist", "aucCoexprDistSameFam" )) {
  
  auc_DT_m <- all_aucCoexprDist_DT[order(all_aucCoexprDist_DT[,coexprdist], decreasing = TRUE),]
  
  if(grepl("old_", coexprdist)) {
    auc_DT_m$datasetLabel <- auc_DT_m$dataset
    auc_DT_m <- auc_DT_m[, c("datasetLabel", "dataset", paste0(coexprdist))]
    auc_DT_m <- unique(auc_DT_m)
  }
  
  auc_DT_m$datasetLabel <- factor(as.character(auc_DT_m$datasetLabel), levels = as.character(auc_DT_m$datasetLabel))
  
  # stopifnot(as.character(auc_DT_m$dataset)  %in% names(dataset_proc_colors) )
  # curr_colors <- dataset_proc_colors[as.character(levels(auc_DT_m$dataset))]
  
  stopifnot(as.character(auc_DT_m$dataset)  %in% names(dataset_proc_colors) )
  curr_colors <- as.character(cancer_subColors[as.character(cancer_subAnnot[as.character(auc_DT_m$dataset)])])
  stopifnot(!is.na(curr_colors))
  
  plotDT <- auc_DT_m
  stopifnot(nrow(plotDT) > 0)
  
  plotDT[,coexprdist] <- 100*(plotDT[,coexprdist] - 1)
  
  colFill <- "deepskyblue4"
  
  myylab <- paste0("% increase CoexprDist AUC")
  myTit <- paste0(var_names[coexprdist])
  my_breaks <- scales::pretty_breaks(n = 10)(plotDT[,coexprdist])
  # my_labels <- my_breaks + 1
  my_labels <- my_breaks
  
  nDS <- length(unique(as.character(plotDT$datasetLabel)))
  mySub <- paste0("(# datasets = ", nDS, ")")
  
  
  p_AUC <- ggplot(plotDT, aes_string(x = paste0("datasetLabel"), y = paste0(coexprdist))) +
    # p_AUC <- ggplot(plotDT, aes_string(x = paste0("dataset"), y = paste0(coexprdist))) +
    geom_bar(stat="identity", position="dodge", width = 0.7, fill = colFill) +
    scale_x_discrete(name="")+
    ggtitle(label=myTit, subtitle = mySub) +
#    scale_y_continuous(name=paste0(myylab),
#                       breaks = my_breaks,
#                       labels = my_labels)+
    scale_y_continuous(name=paste0(myylab),
                       breaks = scales::pretty_breaks(n = 10))+


    coord_cartesian(expand = FALSE) +
#    scale_fill_manual(values = c(auc_fcc = barcolors, auc_coexpr = barcolors),
#                      labels = c(auc_fcc = "FCC", auc_coexpr = "coexpr."),
#                      guide = FALSE )+
    labs(fill  = "") +
    theme( # Increase size of axis lines
      # top, right, bottom and left
      #    plot.margin = unit(c(1, 1, 4.5, 1), "lines"),
      plot.title = element_text(hjust = 0.5, face = "bold", size=16),
      plot.subtitle = element_text(hjust = 0.5),
      panel.grid = element_blank(),
      # panel.grid.major = element_line(colour = "lightpink"),
      # strip.text.x = element_text(size = 6),
      axis.text.x = element_text( hjust=1,vjust = 0.5, size=8, angle = 90, color = curr_colors),
      axis.line.x = element_line(size = .2, color = "black"),
      axis.line.y = element_line(size = .3, color = "black"),
      #    axis.ticks.x = element_blank(),
      axis.text.y = element_text(color="black", hjust=1,vjust = 0.5, size=12),
      axis.title.y = element_text(color="black", size=12),
      axis.title.x = element_text(color="black", size=12),
      panel.border = element_blank(),
      panel.background = element_rect(fill = "transparent"),
      legend.background =  element_rect(),
      legend.key = element_blank()
      # axis.text.x=element_text(size=10, angle=90, vjust=0.5, hjust=1)
    ) #+
  # geom_hline(yintercept = 1, linetype = 2)
  
  if(SSHFS) p_AUC
  
  outFile <- file.path(outFold, paste0(coexprdist, "_TCGA_datasets_barplot.", plotType))
  ggsave(p_AUC, filename = outFile, height = myGGheight, width=myGGwidth)
  cat(paste0("... written: ", outFile, "\n"))
  
}
