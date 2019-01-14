startTime <- Sys.time()

library(foreach)
library(doMC)
library(ggplot2)
library(stringi)

options(scipen=100)

cat("> START: cmp_datasets_MoC.R\n")
# Rscript cmp_datasets_MoC.R

SSHFS <- FALSE
setDir <- ifelse(SSHFS, "/media/electron", "")

source(file.path("utils_fct.R"))
source("../Dixon2018_integrative_data/MoC_heatmap_fct.R")

if(SSHFS) setwd("/media/electron/mnt/etemp/marie/Cancer_HiC_data_TAD_DA")

registerDoMC(ifelse(SSHFS, 2, 40))

outFold <- file.path("CMP_DATASETS_MOC")
system(paste0("mkdir -p ", outFold))

logFile <- file.path(outFold, "cmp_datasets_MoC_logFile.txt")
if(!SSHFS) system(paste0("rm -f ", logFile))
if(SSHFS) logFile <- ""

plotType <- "svg"
widthMoCmat <- 26
heightMoCmat <- 14
widthBoxplot <- 10
heightBoxplot <- 7

strWidthSplit <- 35

binSize <- 40000

txt <- paste0("!! hard-coded bin size:\t", binSize, "\n")
printAndLog(txt, logFile)

consensusPattern <- "_final_domains.txt$" # TO CHECK ZZZ


source("datasets_settings.R")


all_cmps <- combn(all_ds, m = 2)

stopifnot(nrow(all_cmps) == 2)
stopifnot(ncol(all_cmps) > 0)

i=1
all_MoC_dt <- foreach(i = seq_len(ncol(all_cmps)), .combine="rbind") %dopar% {
  
  ds1 <- all_cmps[1,i]
  ds2 <- all_cmps[2,i]

  all_files1 <- eval(parse(text = paste0(ds1, "Files")))
  all_files2 <- eval(parse(text = paste0(ds2, "Files")))
  
  name1 <- eval(parse(text = paste0(ds1, "name")))
  name2 <- eval(parse(text = paste0(ds2, "name")))
  
  cat(paste0("*** START: ", ds1, " vs. ", ds2, "\n"))
    
  chromo = "chr1"
  chr_MoC_dt <- foreach(chromo = intersectChromos, .combine='rbind') %do% {
    
    cat(paste0("> ", ds1, " vs. ", ds2, " - ", chromo, "\n"))
    
    file1 <- all_files1[grepl(paste0(chromo, "_"), basename(all_files1)) & grepl(paste0(name1), all_files1)]
    stopifnot(length(file1) == 1)
    
    file2 <- all_files2[grepl(paste0(chromo, "_"), basename(all_files2)) & grepl(paste0(name2), all_files2)]
    stopifnot(length(file2) == 1)
    
    dt1 <- read.delim(file1, stringsAsFactors = FALSE, header=F, col.names = c("chromo", "start", "end"))
    if(is.character(dt1[1,2]) & is.character(dt1[2,2])) {
      dt1 <- read.delim(file1, stringsAsFactors = FALSE, header=T, col.names = c("chromo", "start", "end"))
    }
    dt2 <- read.delim(file2, stringsAsFactors = FALSE, header=F, col.names = c("chromo", "start", "end"))
    if(is.character(dt2[1,2]) & is.character(dt2[2,2])) {
      dt2 <- read.delim(file2, stringsAsFactors = FALSE, header=T, col.names = c("chromo", "start", "end"))
    }
    stopifnot(ncol(dt1) == 3)
    stopifnot(ncol(dt2) == 3)
    stopifnot(is.numeric(dt1[,2]))
    stopifnot(is.numeric(dt2[,2]))
    stopifnot(is.numeric(dt1[,3]))
    stopifnot(is.numeric(dt2[,3]))
    
    head(dt1, 2)
    head(dt2, 2)
    
    last_before1 <- dt1$end[nrow(dt1)]
    last_after1 <- ceiling(last_before1/binSize)*binSize
    stopifnot(last_after1 >= last_before1)
    dt1$end[nrow(dt1)] <- last_after1
    stopifnot(dt1$end %% binSize == 0)
    stopifnot( (dt1$start-1) %% binSize == 0)
    txt <- paste0("... change last end from\t", last_before1, "\tto\t", last_after1, "\n" )
    printAndLog(txt, logFile)
    
    last_before2 <- dt2$end[nrow(dt2)]
    last_after2 <- ceiling(last_before2/binSize)*binSize
    stopifnot(last_after2 >= last_before2)
    dt2$end[nrow(dt2)] <- last_after2
    stopifnot(dt2$end %% binSize == 0)
    stopifnot( (dt2$start-1) %% binSize == 0)
    txt <- paste0("... change last end from\t", last_before2, "\tto\t", last_after2, "\n" )
    printAndLog(txt, logFile)
    
    stopifnot(!is.na(dt1))
    stopifnot(!is.na(dt2))
    
    chrLen <- max(c(dt1$end, dt2$end))
    
    moc <- calculate_MoC_with_domainTypes(set1=dt1, set2=dt2, chr_len=chrLen,
                                          gaps_as_clusters = TRUE, file_as_input = FALSE)
    
    txt <- paste0("*** MoC\t=\t", round(moc, 4), "\n")
    printAndLog(txt, logFile)
    
    txt <- paste0("\n")
    printAndLog(txt, logFile)
    
    data.frame(
      ds1 = ds1,
      ds2 = ds2,
      chromo = chromo,
      MoC = moc,
      stringsAsFactors = FALSE
    )
  }
  chr_MoC_dt
}

outFile <- file.path(outFold, "all_MoC_dt.Rdata")
save(all_MoC_dt, file = outFile)
cat(paste0("... written: ", outFile, "\n"))

all_MoC_dt <- eval(parse(text = load(outFile)))
stopifnot(nrow(all_MoC_dt) == (length(intersectChromos) * ncol(all_cmps)))


all_MoC_dt_save <- all_MoC_dt

#******************************************************************************************************************************************** DRAW SYMMETRIC MATRIX

# UNCOMMENT TO HAVE NAME WITH CELL LINES ACCESSIONS
#all_MoC_dt$ds1_name <- sapply(all_MoC_dt$ds1, function(x) eval(parse(text = paste0(x, "name"))))
#all_MoC_dt$ds1_init <- all_MoC_dt$ds1
#all_MoC_dt$ds1 <- paste0(all_MoC_dt$ds1, "/\n", all_MoC_dt$ds1_name)

#all_MoC_dt$ds2_name <- sapply(all_MoC_dt$ds2, function(x) eval(parse(text = paste0(x, "name"))))
#all_MoC_dt$ds2_init <- all_MoC_dt$ds2
#all_MoC_dt$ds2 <- paste0(all_MoC_dt$ds2, "/\n", all_MoC_dt$ds2_name)

#all_MoC_dt$ds1 <- gsub("_vs_", " ", all_MoC_dt$ds1)
#all_MoC_dt$ds1 <- unlist(sapply(all_MoC_dt$ds1, function(x) paste0(stri_wrap(str = x, width = strWidthSplit), collapse="\n")))
#all_MoC_dt$ds2 <- gsub("_vs_", " ", all_MoC_dt$ds2)
#all_MoC_dt$ds2 <- unlist(sapply(all_MoC_dt$ds2, function(x) paste0(stri_wrap(str = x, width = strWidthSplit), collapse="\n")))

mean_MoC_dt <- aggregate(MoC ~ ds1 + ds2, FUN=mean, data = all_MoC_dt)
stopifnot(!is.na(mean_MoC_dt))

mean_MoC_dt1 <- mean_MoC_dt[,c("ds1", "ds2", "MoC")]
mean_MoC_dt2 <- mean_MoC_dt[,c("ds2", "ds1", "MoC")]

# self_MoC_dt3 <- data.frame(ds1=all_ds, ds2=all_ds, MoC = 1)
self_MoC_dt3 <- data.frame(ds1=unique(c(all_MoC_dt$ds1, all_MoC_dt$ds2)),
                           ds2=unique(c(all_MoC_dt$ds1, all_MoC_dt$ds2)),
                           MoC = 1)

colnames(mean_MoC_dt2) <- colnames(mean_MoC_dt1)
mocDT <- rbind(rbind( mean_MoC_dt1, mean_MoC_dt2), self_MoC_dt3)
stopifnot(!duplicated(mocDT))

mocDT <- mocDT[order(mocDT$ds1, mocDT$ds2),]
corMat <- reshape(mocDT, idvar="ds1", timevar="ds2", direction="wide")
rownames(corMat) <- corMat$ds1
colnames(corMat) <- sub("MoC.", "", colnames(corMat))
corMat$ds1 <- NULL
corMat <- as.matrix(corMat)
stopifnot(rownames(corMat) == colnames(corMat))
stopifnot(isSymmetric(as.matrix(corMat)))
stopifnot(!is.na(corMat))


tit <- paste0("MoC between tissues (with consensus)\n")

# with pdf output a different unicode character -> save as svg
# outFile <- paste0(outFold, "/", "figure4_MoC_square_heatmap_with_dendro_", curr_norm, "_", res, "kb_all_callers.", "svg")
gplot_dendro <- plot_ggheatmap_with_left_rowdendro(x=as.matrix(corMat),
                                                   ranked_branches =T,
                                                   plotMap = "square", 
                                                   low_limit_col = 0,
                                                   high_limit_col = 1,
                                                   fill_legName = "MoC", 
                                                   dendroLabSize = 4,
                                                   addClusterDot = F,
                                                   annotateMat = TRUE,
                                                   annotateMean = TRUE,
                                                   comparisonName = "caller",
                                                   legCategoryCols = NULL,
                                                   lab_color_vect = NULL)


outFile <- file.path(outFold, paste0("tissues_consensus_MoC_heatmap.", plotType))
ggsave(plot=gplot_dendro, file = outFile, width = widthMoCmat, height = heightMoCmat)
cat(paste0("... written: ", outFile, "\n"))
foo <- try(dev.off())

# outFile <- file.path(outFold, "cmp_MoC_matrix.pdf")
# ggsave(plot=gplot_dendro, file = outFile, width = 26, height = 14)
# cat(paste0("... written: ", outFile, "\n"))

#******************************************************************************************************************************************** BOXPLOT FOR THE CONSENSUS

# all_MoC_dt <- all_MoC_dt_save

# > unique(all_MoC_dt$ds1)

                                                                
# > unique(all_MoC_dt$ds1_init)


cat(paste0(unique(all_MoC_dt$ds1), collapse=","))

for(tissue in c("lung", "breast", "mcf7", "kidney","skin")) {
  
  consensus_dt <- all_MoC_dt[ (grepl(paste0(tissue, "Consensus"), all_MoC_dt$ds1) | grepl(paste0(tissue, "Consensus"), all_MoC_dt$ds2) ) &
                                (grepl(tolower(tissue), tolower(all_MoC_dt$ds1)) & grepl(tolower(tissue), tolower(all_MoC_dt$ds2)) )
                                ,]
  # if short name -> cannot find mcf7 cell lines
  if(nrow(consensus_dt) == 0 & tissue == "mcf7") {
    consensus_dt <- all_MoC_dt[ (grepl(paste0(tissue, "Consensus"), all_MoC_dt$ds1) | grepl(paste0(tissue, "Consensus"), all_MoC_dt$ds2))  &
                                  ( all_MoC_dt$ds1 %in% c("breastCL1", "breastCL2") | all_MoC_dt$ds2 %in% c("breastCL1", "breastCL2") ) 
                                ,]
  }
  # consensus_dt$newDS1 <- ifelse(consensus_dt$ds1 == "consensus", consensus_dt$ds1, consensus_dt$ds2)
  # consensus_dt$newDS2 <- ifelse(consensus_dt$ds2 == "consensus", consensus_dt$ds1, consensus_dt$ds2)
  consensus_dt$newDS1 <- ifelse(grepl(paste0(tissue, "Consensus"), consensus_dt$ds1), consensus_dt$ds1, consensus_dt$ds2)
  consensus_dt$newDS2 <- ifelse(grepl(paste0(tissue, "Consensus"), consensus_dt$ds2), consensus_dt$ds1, consensus_dt$ds2)
        
  stopifnot(grepl(paste0(tissue, "Consensus"), consensus_dt$newDS1))
  stopifnot(!grepl(paste0(tissue, "Consensus"), consensus_dt$newDS2))

  consensus_dt$comp <- paste0(consensus_dt$newDS1, "_", consensus_dt$newDS2)
  stopifnot(sapply(seq_len(nrow(consensus_dt)),function(i) grepl(consensus_dt$ds1[i], consensus_dt$comp[i])))
  stopifnot(sapply(seq_len(nrow(consensus_dt)),function(i) grepl(consensus_dt$ds2[i], consensus_dt$comp[i])))
  
  mean_consensus_dt <- aggregate(MoC ~ newDS1 + newDS2, FUN=mean, data = consensus_dt)
  mean_consensus_dt <- mean_consensus_dt[order(mean_consensus_dt$MoC, decreasing = TRUE),]
  consensus_dt$newDS2 <- factor(as.character(consensus_dt$newDS2), levels = mean_consensus_dt$newDS2)
  consensus_dt$chromo <- factor(as.character(consensus_dt$chromo), levels = paste0("chr", c(1:22, "X")))
  
  stopifnot(!is.na(consensus_dt))
  
  p_common <- ggplot(consensus_dt, aes(x = newDS2, y = MoC)) + 
    geom_boxplot(outlier.shape=NA) +
    # geom_jitter(aes(colour = chromo)) +
    scale_x_discrete(name="")+
    # scale_y_continuous(name=paste0("-log10(", padjVarGO, ")"),
    scale_y_continuous(name=paste0("MoC with consensus"),
                       breaks = scales::pretty_breaks(n = 10))+ #, limits = c(0, max(auc_DT_m$value)+0.05))+
    # coord_cartesian(expand = FALSE) +
    # scale_fill_manual(values = c(selectGenes = "dodgerblue4", selectTADs_genes = "darkorange2"),
    #                   labels = c(selectGenes = "selectGenes", selectTADs_genes = "selectTADs_genes"))+
    # scale_colour_manual(values = c(selectGenes = "dodgerblue4", selectTADs_genes = "darkorange2"),
    #                     labels = c(selectGenes = "selectGenes", selectTADs_genes = "selectTADs_genes"), guide = F)+
    labs(colour  = "") +
    ggtitle(label = paste0(tissue, " - MoC with consensus"))+
    theme( # Increase size of axis lines
      # top, right, bottom and left
      # plot.margin = unit(c(1, 1, 4.5, 1), "lines"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size=10),
      panel.grid = element_blank(),
      # panel.grid.major = element_line(colour = "lightpink"),
      # strip.text.x = element_text(),
      axis.text.x = element_text( hjust=1,vjust = 0.5, size=12, angle = 90),
      axis.line.x = element_line(size = .2, color = "black"),
      axis.line.y = element_line(size = .3, color = "black"),
      #    axis.ticks.x = element_blank(),
      axis.text.y = element_text(color="black", hjust=1,vjust = 0.5),
      axis.title.y = element_text(color="black", size=12),
      axis.title.x = element_text(color="black", size=12),
      panel.border = element_blank(),
      panel.background = element_rect(fill = "transparent"),
      legend.background =  element_rect(),
      legend.key = element_blank()
      # axis.text.x=element_text(size=10, angle=90, vjust=0.5, hjust=1)
    ) #+
  # geom_hline(yintercept = 1, linetype = 2)
  
  if(SSHFS) p_all
  
  p_dot <- p_common +  geom_jitter(aes(colour = chromo))
  if(SSHFS) p_dot
  
  p_txt <- p_common + geom_text(aes(label=chromo, colour=chromo, fontface="bold"),size=2.5, position = position_jitter(w = 0.3)) + guides(colour = "none")
  if(SSHFS) p_txt
  
  
  outFile <- file.path(outFold, paste0("MoC_", tissue, "_consensus_boxplot_chromoDots.", plotType))
  ggsave(plot=p_dot, file = outFile, width = widthBoxplot, height = heightBoxplot)
  cat(paste0("... written: ", outFile, "\n"))
  foo <- try(dev.off())
  
  outFile <- file.path(outFold, paste0("MoC_matching_", tissue, "_consensus_boxplot_chromoLabs.", plotType))
  ggsave(plot=p_txt, file = outFile, width = widthBoxplot, height = heightBoxplot)
  cat(paste0("... written: ", outFile, "\n"))
  foo <- try(dev.off())
  
  
}


### add the same boxplot here for pipTopDomconsensus only !!!

tissue <- pipConsensusname
if( all(all_MoC_dt$ds1 %in% all_ds)) tissue <- "pipConsensus"

consensus_dt <- all_MoC_dt[ (grepl(paste0(tissue), all_MoC_dt$ds1) | grepl(paste0(tissue), all_MoC_dt$ds2) ),]
# consensus_dt$newDS1 <- ifelse(consensus_dt$ds1 == "consensus", consensus_dt$ds1, consensus_dt$ds2)
# consensus_dt$newDS2 <- ifelse(consensus_dt$ds2 == "consensus", consensus_dt$ds1, consensus_dt$ds2)
consensus_dt$newDS1 <- ifelse(grepl(paste0(tissue), consensus_dt$ds1), consensus_dt$ds1, consensus_dt$ds2)
consensus_dt$newDS2 <- ifelse(grepl(paste0(tissue), consensus_dt$ds2), consensus_dt$ds1, consensus_dt$ds2)
    
stopifnot(grepl(paste0(tissue), consensus_dt$newDS1))
stopifnot(!grepl(paste0(tissue), consensus_dt$newDS2))

consensus_dt$comp <- paste0(consensus_dt$newDS1, "_", consensus_dt$newDS2)
stopifnot(sapply(seq_len(nrow(consensus_dt)),function(i) grepl(consensus_dt$ds1[i], consensus_dt$comp[i])))
stopifnot(sapply(seq_len(nrow(consensus_dt)),function(i) grepl(consensus_dt$ds2[i], consensus_dt$comp[i])))

mean_consensus_dt <- aggregate(MoC ~ newDS1 + newDS2, FUN=mean, data = consensus_dt)
mean_consensus_dt <- mean_consensus_dt[order(mean_consensus_dt$MoC, decreasing = TRUE),]
consensus_dt$newDS2 <- factor(as.character(consensus_dt$newDS2), levels = mean_consensus_dt$newDS2)
consensus_dt$chromo <- factor(as.character(consensus_dt$chromo), levels = paste0("chr", c(1:22, "X")))

stopifnot(!is.na(consensus_dt))

p_common <- ggplot(consensus_dt, aes(x = newDS2, y = MoC)) + 
geom_boxplot(outlier.shape=NA) +
# geom_jitter(aes(colour = chromo)) +
scale_x_discrete(name="")+
# scale_y_continuous(name=paste0("-log10(", padjVarGO, ")"),
scale_y_continuous(name=paste0("MoC"),
                   breaks = scales::pretty_breaks(n = 10))+ #, limits = c(0, max(auc_DT_m$value)+0.05))+
# coord_cartesian(expand = FALSE) +
# scale_fill_manual(values = c(selectGenes = "dodgerblue4", selectTADs_genes = "darkorange2"),
#                   labels = c(selectGenes = "selectGenes", selectTADs_genes = "selectTADs_genes"))+
# scale_colour_manual(values = c(selectGenes = "dodgerblue4", selectTADs_genes = "darkorange2"),
#                     labels = c(selectGenes = "selectGenes", selectTADs_genes = "selectTADs_genes"), guide = F)+
labs(colour  = "") +
ggtitle(label = paste0("MoC with ", tissue))+
theme( # Increase size of axis lines
  # top, right, bottom and left
  # plot.margin = unit(c(1, 1, 4.5, 1), "lines"),
  plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
  plot.subtitle = element_text(hjust = 0.5, size=10),
  panel.grid = element_blank(),
  # panel.grid.major = element_line(colour = "lightpink"),
  # strip.text.x = element_text(),
  axis.text.x = element_text( hjust=1,vjust = 0.5, size=12, angle = 90),
  axis.line.x = element_line(size = .2, color = "black"),
  axis.line.y = element_line(size = .3, color = "black"),
  #    axis.ticks.x = element_blank(),
  axis.text.y = element_text(color="black", hjust=1,vjust = 0.5),
  axis.title.y = element_text(color="black", size=12),
  axis.title.x = element_text(color="black", size=12),
  panel.border = element_blank(),
  panel.background = element_rect(fill = "transparent"),
  legend.background =  element_rect(),
  legend.key = element_blank()
  # axis.text.x=element_text(size=10, angle=90, vjust=0.5, hjust=1)
) #+
# geom_hline(yintercept = 1, linetype = 2)

if(SSHFS) p_all

p_dot <- p_common +  geom_jitter(aes(colour = chromo))
if(SSHFS) p_dot

p_txt <- p_common + geom_text(aes(label=chromo, colour=chromo, fontface="bold"),size=2.5, position = position_jitter(w = 0.3)) + guides(colour = "none")
if(SSHFS) p_txt


outFile <- file.path(outFold, paste0("MoC_matching_", tissue, "_consensus_boxplot_chromoDots.", plotType))
ggsave(plot=p_dot, file = outFile, width = widthBoxplot, height = heightBoxplot)
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFold, paste0("MoC_matching_", tissue, "_consensus_boxplot_chromoLabs.", plotType))
ggsave(plot=p_txt, file = outFile, width = widthBoxplot, height = heightBoxplot)
cat(paste0("... written: ", outFile, "\n"))


######################################################################################
cat(paste0("... written: ", logFile, "\n"))
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))






