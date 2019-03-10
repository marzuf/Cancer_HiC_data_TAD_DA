startTime <- Sys.time()

library(foreach)
library(doMC)
library(ggplot2)

options(scipen=100)

# Rscript cmp_auc_fcc_tissue.R

cat("> START: cmp_auc_fcc_tissue.R\n")

SSHFS <- FALSE
# setDir <- ifelse(SSHFS, "/media/electron", "")
setDir <- ifelse(SSHFS, "~/media/electron", "")


script17name <- "170revision2EZH2_score_auc_pval_permGenes"

plotType <- "svg"
myHeight <- ifelse(plotType == "png", 300, 7)
myWidth <- myHeight
myGGheight <- ifelse(plotType == "png", 300, 7)
myGGwidth <- ifelse(plotType == "png", 500, 10)
myGGwidth <- ifelse(plotType == "png", 500, 10)
plotCex <- 1.2

outFold <- file.path("CMP_AUC_FCC_TISSUE")
system(paste("mkdir -p ", outFold))

source( file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg"), "set_dataset_colors.R"))
head(score_DT)

source( file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_11_18"), "analysis_utils.R"))
source( file.path("colors_utils.R"))


dataset_proc_colors <- setNames(score_DT$proc_col, score_DT$dataset)
length(dataset_proc_colors)

my_colors_leg <- my_colors

#### BUILD THE BARPLOT AUC FCC RANKING 

# retrieve all the FCC auc ratios available


all_ratio_files <- list.files(file.path("PIPELINE", "OUTPUT_FOLDER"),
                          pattern = "auc_ratios.Rdata",
                          recursive=TRUE,
                          full.names=TRUE)
stopifnot(length(all_ratio_files) > 0)

curr_file = file.path(setDir,
                      "/mnt/etemp/marie/Cancer_HiC_data_TAD_DA",
                      "PIPELINE/OUTPUT_FOLDER/ENCSR079VIJ_G401_40kb",
                        "TCGAkich_norm_kich",
                        script17name, "auc_ratios.Rdata")

all_aucFCC_DT <- foreach(curr_file = all_ratio_files, .combine='rbind') %dopar% {
  curr_ds <- basename(dirname(dirname(curr_file)))
  stopifnot(grepl("TCGA", curr_ds))
  all_ratios <- eval(parse(text = load(curr_file)))
  aucfcc <- as.numeric(all_ratios["prodSignedRatio_auc_permGenes"])
  stopifnot(! is.na(aucfcc))
  
  
  data.frame(
    datasetLabel = paste0(basename(dirname(dirname(dirname(curr_file)))), "\n", basename(dirname(dirname(curr_file)))),
    dataset = curr_ds,
    aucFCC = aucfcc,
    stringsAsFactors = FALSE
  )
}

myTit <- paste0(" % increase FCC AUC\nTissue-specific consensus vs. pipeline consensus")
#myTit <- paste0("Tissue-specific consensus vs. pipeline consensus % increase FCC AUC")
x_lab <- paste0("% increase FCC AUC - pipeline consensus (pipeline consensus)")
y_lab <- paste0("% increase FCC AUC - tissue consensus (tissue consensus)")
mySub <- paste0("(# DS = ", nrow(all_aucFCC_DT), ")")

myy <- 100*(all_aucFCC_DT$aucFCC-1)
mynames <- all_aucFCC_DT$dataset


curr_colors <- as.character(cancer_subColors[as.character(cancer_subAnnot[mynames])])
stopifnot(!is.na(curr_colors))




##############################################################################################
############################################################################################## BARPLOT WITH ONE AUC RATIOS
##############################################################################################

var_names <- c(aucFCC = "% increase AUC FCC (tissue-specific)",
               old_aucFCC = "% increase AUC FCC (pipeline consensus)"
               )

for(fcc in c("aucFCC")) {
  
  
  auc_DT_m <- all_aucFCC_DT[order(all_aucFCC_DT[,fcc], decreasing = TRUE),]
  
  if(fcc == "old_aucFCC") {
    auc_DT_m$datasetLabel <- auc_DT_m$dataset
    auc_DT_m <- auc_DT_m[, c("dataset", "datasetLabel", paste0(fcc))]
    auc_DT_m <- unique(auc_DT_m)
  }
  
  # auc_DT_m$dataset <- factor(as.character(auc_DT_m$dataset), levels = unique(as.character(auc_DT_m$dataset)))
  auc_DT_m$datasetLabel <- factor(as.character(auc_DT_m$datasetLabel), levels = as.character(auc_DT_m$datasetLabel))
  
  
  # stopifnot(as.character(auc_DT_m$dataset)  %in% names(dataset_proc_colors) )
  # curr_colors <- dataset_proc_colors[as.character(levels(auc_DT_m$dataset))]
  
  cat(paste0(as.character(auc_DT_m$dataset)[!as.character(auc_DT_m$dataset)  %in% names(dataset_proc_colors) ], collapse="\n"), "\n")
  
  stopifnot(as.character(auc_DT_m$dataset)  %in% names(dataset_proc_colors) )
  
  curr_colors <- as.character(cancer_subColors[as.character(cancer_subAnnot[as.character(auc_DT_m$dataset)])])
  stopifnot(!is.na(curr_colors))
  
  plotDT <- auc_DT_m
  stopifnot(nrow(plotDT) > 0)
  
  plotDT[,fcc] <- 100*(plotDT[,fcc] - 1)
  
  colFill <- "deepskyblue4"
  
  myylab <- paste0("% FCC AUC increase")
  myTit <- paste0(var_names[fcc])
  my_breaks <- scales::pretty_breaks(n = 10)(plotDT[,fcc])
  # my_labels <- my_breaks + 1
  my_labels <- my_breaks
  
  nDS <- length(unique(as.character(plotDT$datasetLabel)))
  mySub <- paste0("(# datasets = ", nDS, ")")
  
  # p_AUC <- ggplot(plotDT, aes_string(x = paste0("dataset"), y = paste0(fcc))) +
    p_AUC <- ggplot(plotDT, aes_string(x = paste0("datasetLabel"), y = paste0(fcc))) +
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
  
  outFile <- file.path(outFold, paste0(fcc, "_TCGA_datasets_barplot.", plotType))
  ggsave(p_AUC, filename = outFile, height = myGGheight, width=myGGwidth)
  cat(paste0("... written: ", outFile, "\n"))
  
}
