startTime <- Sys.time()

library(foreach)
library(doMC)
library(ggplot2)
library(stringi)

options(scipen=100)

cat("> START: cmp_auc_coexprDist_pipeline_tissue.R\n")
# Rscript cmp_auc_coexprDist_pipeline_tissue.R

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


plotType <- "png"
myHeight <- ifelse(plotType == "png", 300, 7)
myWidth <- myHeight
myGGheight <- ifelse(plotType == "png", 300, 7)
myGGwidth <- ifelse(plotType == "png", 500, 10)
plotCex <- 1.2

outFold <- file.path("CMP_AUC_COEXPRDIST_PIPELINE_TISSUE")
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

all_aucCoexprDist_DT <- foreach(curr_file = all_ratio_files, .combine='rbind') %dopar% {
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
  
  data.frame(
    dataset = ds_name,
    aucCoexprDist = aucCoexprDist,
    old_aucCoexprDist = old_aucCoexprDist,
    aucCoexprDistSameFam = aucCoexprDistSameFam,
    old_aucCoexprDistSameFam = old_aucCoexprDistSameFam, 
    stringsAsFactors = FALSE
  )
}


#################################### 
### correlation old and new scores  -- COEXPR DIST
#################################### 

myTit <- paste0(" % increase CoexprDist AUC\nTissue-specific consensus vs. pipeline consensus")
#myTit <- paste0("Tissue-specific consensus vs. pipeline consensus % increase FCC AUC")
x_lab <- paste0("% increase CoexprDist AUC - pipeline consensus (pipeline consensus)")
y_lab <- paste0("% increase CoexprDist AUC - tissue consensus (tissue consensus)")
mySub <- paste0("(# DS = ", nrow(all_aucCoexprDist_DT), ")")

myx <- 100*(all_aucCoexprDist_DT$old_aucCoexprDist-1)
myy <- 100*(all_aucCoexprDist_DT$aucCoexprDist-1)
mynames <- all_aucCoexprDist_DT$dataset


curr_colors <- as.character(cancer_subColors[as.character(cancer_subAnnot[mynames])])
stopifnot(!is.na(curr_colors))

outFile <- file.path(outFold, paste0("pipeline_vs_tissue_coexprDistIncrease.", plotType))
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
legend("topleft",
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


#################################### 
### correlation old and new scores  -- COEXPR DIST SAME FAM
#################################### 

myTit <- paste0(" % increase CoexprDist AUC (same fam.)\nTissue-specific consensus vs. pipeline consensus")
#myTit <- paste0("Tissue-specific consensus vs. pipeline consensus % increase FCC AUC")
x_lab <- paste0("% increase CoexprDist AUC - pipeline consensus (pipeline consensus)")
y_lab <- paste0("% increase CoexprDist AUC - tissue consensus (tissue consensus)")
mySub <- paste0("(# DS = ", nrow(all_aucCoexprDist_DT), ")")

myx <- 100*(all_aucCoexprDist_DT$old_aucCoexprDistSameFam-1)
myy <- 100*(all_aucCoexprDist_DT$aucCoexprDistSameFam-1)
mynames <- all_aucCoexprDist_DT$dataset

curr_colors <- as.character(cancer_subColors[as.character(cancer_subAnnot[mynames])])
stopifnot(!is.na(curr_colors))

outFile <- file.path(outFold, paste0("pipeline_vs_tissue_coexprDistSameFamIncrease.", plotType))
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
legend("topleft",
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


##############################################################################################
############################################################################################## BARPLOT WITH ONE AUC RATIOS
##############################################################################################

var_names <- c(aucCoexprDist = "% increase AUC coexprDist (tissue specific)",
               old_aucCoexprDist = "% increase AUC coexprDist (pipeline consensus)",
               aucCoexprDistSameFam = "% increase AUC coexprDist (same fam.; tissue specific)",
               old_aucCoexprDistSameFam = "% increase AUC coexprDist (same fam.; pipeline consensus)"
               )

for(coexprdist in c("aucCoexprDist", "old_aucCoexprDist", "aucCoexprDistSameFam", "old_aucCoexprDistSameFam" )) {
  
  auc_DT_m <- all_aucCoexprDist_DT[order(all_aucCoexprDist_DT[,coexprdist], decreasing = TRUE),]
  
  auc_DT_m$dataset <- factor(as.character(auc_DT_m$dataset), levels = as.character(auc_DT_m$dataset))
  
  # stopifnot(as.character(auc_DT_m$dataset)  %in% names(dataset_proc_colors) )
  # curr_colors <- dataset_proc_colors[as.character(levels(auc_DT_m$dataset))]
  
  stopifnot(as.character(auc_DT_m$dataset)  %in% names(dataset_proc_colors) )
  curr_colors <- as.character(cancer_subColors[as.character(cancer_subAnnot[levels(auc_DT_m$dataset)])])
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
  
  p_AUC <- ggplot(plotDT, aes_string(x = paste0("dataset"), y = paste0(coexprdist))) +
    geom_bar(stat="identity", position="dodge", width = 0.7, fill = colFill) +
    scale_x_discrete(name="")+
    ggtitle(label=myTit) +
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
      panel.grid = element_blank(),
      # panel.grid.major = element_line(colour = "lightpink"),
      # strip.text.x = element_text(size = 6),
      axis.text.x = element_text( hjust=1,vjust = 0.5, size=10, angle = 90, color = curr_colors),
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
