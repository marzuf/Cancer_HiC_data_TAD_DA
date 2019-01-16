startTime <- Sys.time()

library(foreach)
library(doMC)
library(ggplot2)

options(scipen=100)

# !!! need to have run check_matResol.R before !!!

cat("> START: cmp_datasets_resol.R\n")
# Rscript cmp_datasets_resol.R

SSHFS <- FALSE
setDir <- ifelse(SSHFS, "/media/electron", "")
setDir <- ifelse(SSHFS, "~/media/electron", "")

source(file.path("utils_fct.R"))
source("datasets_settings.R")

# if(SSHFS) setwd("/media/electron/mnt/etemp/marie/Cancer_HiC_data_TAD_DA")
if(SSHFS) setwd("~/media/electron/mnt/etemp/marie/Cancer_HiC_data_TAD_DA")


registerDoMC(ifelse(SSHFS, 2, 40))

outFold <- file.path("CMP_DATASETS_RESOL")
dir.create(outFold, recursive = TRUE)

logFile <- file.path(outFold, "cmp_datasets_resol_logFile.txt")
if(!SSHFS) file.remove(logFile)
if(SSHFS) logFile <- ""

plotType <- "svg"
plotType <- "png"

widthBoxplot <- ifelse(plotType =="png", 600, 10)
heightBoxplot <- ifelse(plotType =="png", 400, 7)

widthBoxplot <- 10
heightBoxplot <- 7

binSize <- 40000
binSizeKb <- binSize/1000

all_chromo <- paste0("chr", c(1:22, "X"))  # for plot ordering

txt <- paste0("!! hard-coded bin size:\t", binSize, "\n")
printAndLog(txt, logFile)

mainFold <- "CHECK_MATRESOL"
stopifnot(dir.exists(mainFold))

# retrieve the intersect chromos:
all_files_rowSum <- list.files(mainFold, full.names=T, recursive = TRUE, pattern="_matrixRowSum\\.Rdata$")
stopifnot(length(all_files_rowSum) > 0)
files_chr_dt <- data.frame(
  ds = basename(dirname(all_files_rowSum)),
  chromo = as.character(sapply(basename(all_files_rowSum), function(x) gsub(".+_(chr.+?)_matrixRowSum.Rdata", "\\1", x))),
  stringsAsFactors = FALSE
)
# take intersect chromos here
intersectChromos <- Reduce(intersect, by(files_chr_dt$ds, data=files_chr_dt, FUN=function(x) unique(as.character(x$chromo))))
txt <- paste0("... # intersectChromos:\t", length(intersectChromos), "\n")
printAndLog(txt, logFile)
txt <- paste0("... intersectChromos:\t", paste0(intersectChromos, collapse=","), "\n")
printAndLog(txt, logFile)

# for each  chromo -> row sum boxplot
  var_names <- c(
    "rowSum"= "Hi-C matrix row sum", 
    "rowSum_log10"= "Hi-C matrix row sum [log10]", 
    "rowSumNoOut"= "Hi-C matrix row sum (5-95% quantile)", 
    "rowSumNoOut_log10" = "Hi-C matrix row sum (5-95% quantile) [log10]", 
     "matrixDim"= "Hi-C matrix dimension"
  )

chromo="chr1"
all_chromo_DT <- foreach(chromo = intersectChromos, .combine="rbind") %dopar% {
  
  chromo_files <- all_files_rowSum[grepl(paste0("_", chromo, "_"), all_files_rowSum)]
  
  chromo_DT <- foreach(chr_file = chromo_files, .combine='rbind') %do% {
    
    curr_ds <- gsub("^(.+)_(chr.+)_matrixRowSum.Rdata$", "\\1", basename(chr_file))
    
    stopifnot(length(curr_ds) > 0)
    
    rowsum_chromo <- eval(parse(text = load(chr_file)))
    
    lowOut <- as.numeric(quantile(rowsum_chromo, probs=0.05))
    highOut <- as.numeric(quantile(rowsum_chromo, probs=0.95))
    
    rowsum_chromo_noout <- rowsum_chromo
    rowsum_chromo_noout[rowsum_chromo_noout <= lowOut | rowsum_chromo_noout >= highOut] <- NA
    
    stopifnot(length(rowsum_chromo) > 0)
    
    stopifnot(length(rowsum_chromo) == length(rowsum_chromo_noout) )
    
    matrixDim <- length(rowsum_chromo)
    
    data.frame(
      dataset=curr_ds,
      chromo = chromo,
      matrixDim = matrixDim,
      rowIdx = 1:length(rowsum_chromo),
      rowSum = rowsum_chromo,
      rowSum_log10 = log10(rowsum_chromo),
      rowSumNoOut = rowsum_chromo_noout,
      rowSumNoOut_log10 = log10(rowsum_chromo_noout),
      stringsAsFactors = FALSE
    )
  }# end iterate over datasets
  
  stopifnot(chromo_DT$dataset %in% cl_names)
  
  chromo_DT$dataset_label <- unlist(sapply(as.character(chromo_DT$dataset), function(x) {
    dslab <- as.character(cl_labs[x])
    if(is.na(dslab)) {
      paste0(names(cl_names[cl_names == x]))
    } else {
      paste0(names(cl_names[cl_names == x]), "\n(", dslab, ")")
    }
  } ))

  for(var_to_plot in c("rowSum", "rowSum_log10", "rowSumNoOut",  "rowSumNoOut_log10", "matrixDim")) {
    
    var_tit <- var_names[var_to_plot]
    stopifnot(!is.na(var_tit))
    
    tmp_DT <- aggregate(as.formula(paste0(var_to_plot, " ~ dataset_label")), data = chromo_DT, FUN=mean, na.rm=T)
    tmp_DT <- tmp_DT[order(tmp_DT[, var_to_plot], decreasing=TRUE),]
    
    plot_chromo_DT <- chromo_DT[!is.infinite(chromo_DT[,var_to_plot]),]
    plot_chromo_DT$dataset_label <- factor(as.character(plot_chromo_DT$dataset_label), levels = as.character(tmp_DT$dataset_label))
    
    p_chromo_all_ds <- ggplot(plot_chromo_DT, aes_string(x = "dataset_label", y = var_to_plot))+
      geom_boxplot()+
      ggtitle(paste0(chromo, ": ", var_tit)) + 
      scale_x_discrete(name="")+
      scale_y_continuous(name = paste0(var_tit), breaks = scales::pretty_breaks(n = 10))+ #, limits = c(0, max(auc_DT_m$value)+0.05))+
      labs(colour  = "") +
      theme( # Increase size of axis lines
        # top, right, bottom and left
        # plot.margin = unit(c(1, 1, 4.5, 1), "lines"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, size=10),
        panel.grid = element_blank(),
        axis.text.x = element_text( hjust=1,vjust = 0.5, size=14, angle = 90),
        axis.line.x = element_line(size = .2, color = "black"),
        axis.line.y = element_line(size = .3, color = "black"),
        axis.text.y = element_text(color="black", hjust=1,vjust = 0.5),
        axis.title.y = element_text(color="black", size=12),
        axis.title.x = element_text(color="black", size=12),
        panel.border = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        legend.background =  element_rect(),
        legend.key = element_blank()
      )
    outFile <- file.path(outFold, paste0(var_to_plot, "_", chromo, "_boxplot.", plotType))
    ggsave(plot=p_chromo_all_ds, file = outFile, width = widthBoxplot, height = heightBoxplot)
    cat(paste0("... written: ", outFile, "\n"))
  }
  # stopifnot(!is.na(chromo_DT))
  chromo_DT
} # end iterate over chromos

outFile <- file.path(outFold, "all_chromo_DT.Rdata")
save(all_chromo_DT, file = outFile)
cat(paste0("... written: ", outFile, "\n"))
# stopifnot(!is.na(all_chromo_DT))

my_ggtit <- "all chromos"
var_to_plot <- "rowSum_log10"

for(var_to_plot in c("rowSum", "rowSum_log10", "rowSumNoOut", "rowSumNoOut_log10", "matrixDim")) {
  
  var_tit <- var_names[var_to_plot]
  stopifnot(!is.na(var_tit))
  
  tmp_DT <- all_chromo_DT[!is.infinite(all_chromo_DT[,var_to_plot]),]
  tmp_DT <- aggregate(as.formula(paste0(var_to_plot, " ~ dataset_label")), data = tmp_DT, FUN=mean, na.rm=T)
  tmp_DT <- tmp_DT[order(tmp_DT[, var_to_plot], decreasing=TRUE),]
  
  plot_all_chromo_DT <- all_chromo_DT[!is.infinite(all_chromo_DT[,var_to_plot]),]
  plot_all_chromo_DT$dataset_label <- factor(as.character(plot_all_chromo_DT$dataset_label), levels = as.character(tmp_DT$dataset_label))
  
  p_all_chromo_all_ds <- ggplot(plot_all_chromo_DT, aes_string(x = "dataset_label", y = var_to_plot))+
    geom_boxplot()+
    # ggtitle(paste0(var_to_plot, ": ", my_ggtit)) +
    ggtitle(paste0(var_tit, ": ", my_ggtit)) +
    scale_x_discrete(name="")+
    scale_y_continuous(name=paste0(var_tit), breaks = scales::pretty_breaks(n = 10))+
    labs(colour  = "") +
    theme( # Increase size of axis lines
      # top, right, bottom and left
      # plot.margin = unit(c(1, 1, 4.5, 1), "lines"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size=10),
      panel.grid = element_blank(),
      axis.text.x = element_text( hjust=1,vjust = 0.5, size=14, angle = 90),
      axis.line.x = element_line(size = .2, color = "black"),
      axis.line.y = element_line(size = .3, color = "black"),
      axis.text.y = element_text(color="black", hjust=1,vjust = 0.5),
      axis.title.y = element_text(color="black", size=12),
      axis.title.x = element_text(color="black", size=12),
      panel.border = element_blank(),
      panel.background = element_rect(fill = "transparent"),
      legend.background =  element_rect(),
      legend.key = element_blank()
    )
  outFile <- file.path(outFold, paste0("all_chromos_all_datasets_boxplot_", var_to_plot, ".", plotType))
  ggsave(plot=p_all_chromo_all_ds, file = outFile, width = widthBoxplot, height = heightBoxplot)
  cat(paste0("... written: ", outFile, "\n"))
}

all_files_resol <- list.files(mainFold, full.names=T, recursive = TRUE, pattern="_check_resolDT\\.Rdata$")
stopifnot(length(all_files_resol) > 0)

curr_file <- all_files_resol[1]
all_resol_DT <- foreach(curr_file = all_files_resol, .combine='rbind') %do% {
  eval(parse(text=load(curr_file)))
}
all_resol_DT$countSum_log10 <- log10(all_resol_DT$countSum)
all_vars <- colnames(all_resol_DT)[!colnames(all_resol_DT) %in% c("dataset", "chromo")]
all_resol_DT$datasetLabel <- unlist(sapply(as.character(all_resol_DT$dataset), function(x) names(cl_names[cl_names == x])))

var_names["countSum"] <- "Hi-C matrix tot. contact sum"
var_names["countSum_log10"] <- "Hi-C matrix tot. contact sum [log10]"
var_names["rowAbove1000"] <- "Ratio of rows with > 1000 counts"

outFile <- file.path(outFold, "all_resol_DT.Rdata")
save(all_resol_DT, file = outFile)
cat(paste0("... written: ", outFile, "\n"))
stopifnot(!is.na(all_resol_DT))

var_to_plot = "rowAbove1000"
for(var_to_plot in all_vars) {
  
  var_tit <- var_names[var_to_plot]
  stopifnot(!is.na(var_tit))
  
  mean_all_resol_DT <- aggregate(as.formula(paste0(var_to_plot, " ~ dataset")), FUN=mean, data = all_resol_DT)
  mean_all_resol_DT <- mean_all_resol_DT[order(mean_all_resol_DT[, var_to_plot], decreasing = TRUE),]
  all_resol_DT$dataset <- factor(as.character(all_resol_DT$dataset), levels = mean_all_resol_DT$dataset)
  
  mean_all_resol_DT <- aggregate(as.formula(paste0(var_to_plot, " ~ datasetLabel")), FUN=mean, data = all_resol_DT)
  mean_all_resol_DT <- mean_all_resol_DT[order(mean_all_resol_DT[, var_to_plot], decreasing = TRUE),]
  all_resol_DT$datasetLabel <- factor(as.character(all_resol_DT$datasetLabel), levels = mean_all_resol_DT$datasetLabel)
  
  # p_common <- ggplot(all_resol_DT, aes_string(x = "dataset", y = var_to_plot)) + 
  p_common <- ggplot(all_resol_DT, aes_string(x = "datasetLabel", y = var_to_plot)) + 
    geom_boxplot(outlier.shape=NA) +
    scale_x_discrete(name="")+
    scale_y_continuous(name=paste0(var_tit),
                       breaks = scales::pretty_breaks(n = 10))+ 
    labs(colour  = "") +
    ggtitle(label = paste0(var_tit))+
    theme( # Increase size of axis lines
      # top, right, bottom and left
      # plot.margin = unit(c(1, 1, 4.5, 1), "lines"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size=10),
      panel.grid = element_blank(),
      axis.text.x = element_text( hjust=1,vjust = 0.5, size=14, angle = 90),
      axis.line.x = element_line(size = .2, color = "black"),
      axis.line.y = element_line(size = .3, color = "black"),
      axis.text.y = element_text(color="black", hjust=1,vjust = 0.5),
      axis.title.y = element_text(color="black", size=12),
      axis.title.x = element_text(color="black", size=12),
      panel.border = element_blank(),
      panel.background = element_rect(fill = "transparent"),
      legend.background =  element_rect(),
      legend.key = element_blank()
    )
  
  if(SSHFS) p_all
  p_dot <- p_common +  geom_jitter(aes(colour = chromo))
  if(SSHFS) p_dot
  p_txt <- p_common + geom_text(aes(label=chromo, colour=chromo, fontface="bold"),size=2.5, position = position_jitter(w = 0.3)) + guides(colour = "none")
  if(SSHFS) p_txt
  # outFile <- file.path(outFold, paste0(var_to_plot, "_boxplot_chromoDots.", plotType))
  # ggsave(plot=p_dot, file = outFile, width = widthBoxplot, height = heightBoxplot)
  # cat(paste0("... written: ", outFile, "\n"))
  outFile <- file.path(outFold, paste0(var_to_plot, "_boxplot_chromoLabs.", plotType))
  ggsave(plot=p_txt, file = outFile, width = widthBoxplot, height = heightBoxplot)
  cat(paste0("... written: ", outFile, "\n"))
}



######################################################################################
cat(paste0("... written: ", logFile, "\n"))
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))






