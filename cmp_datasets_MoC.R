startTime <- Sys.time()

library(foreach)
library(doMC)
library(ggplot2)

options(scipen=100)

cat("> START: cmp_datasets_MoC.R\n")
# Rscript cmp_datasets_MoC.R

buildTable <- TRUE

SSHFS <- FALSE
setDir <- ifelse(SSHFS, "/media/electron", "")
setDir <- ifelse(SSHFS, "~/media/electron", "")

# if(SSHFS) setwd("/media/electron/mnt/etemp/marie/Cancer_HiC_data_TAD_DA")
if(SSHFS) setwd("~/media/electron/mnt/etemp/marie/Cancer_HiC_data_TAD_DA")

source(file.path("utils_fct.R"))
source("../Dixon2018_integrative_data/MoC_heatmap_fct.R")
source(file.path("datasets_settings.R"))

registerDoMC(ifelse(SSHFS, 2, 40))

outFold <- file.path("CMP_DATASETS_MOC")
dir.create(outFold)

logFile <- file.path(outFold, "cmp_datasets_MoC_logFile.txt")
if(!SSHFS) file.remove(logFile)
if(SSHFS) logFile <- ""

plotType <- "png"
widthMoCmat <- 26
heightMoCmat <- 14
widthBoxplot <- 10
heightBoxplot <- 7

binSize <- 40000
binSizeKb <- binSize/1000

txt <- paste0("!! hard-coded buildTable:\t", as.character(buildTable), "\n")
printAndLog(txt, logFile)

txt <- paste0("!! hard-coded bin size:\t", binSize, "\n")
printAndLog(txt, logFile)

folderSuffix <- paste0("_", binSizeKb, "kb")
TADfilePattern <- "_final_domains.txt$"
TADfolder <- "FINAL_DOMAINS"

all_chromos <- c(paste0("chr", 1:23), "chrX")  # also used for ordering on the plots

xlabType <- "_label"
stopifnot(xlabType %in% c("", "_label"))

### SELECT WHICH CELL LINES TO INCLUDE IN THE COMPARISONS
cl_to_cmp <- c(
  "MCF-7",
  "ENCSR549MGQ_T47D",
  "MCF-7ENCSR549MGQ_T47D",
  
  "GSE105318_DLD1",                              # !!! RUNNING !!!
  
  "ENCSR079VIJ_G401",
  "ENCSR401TBQ_Caki2",
  "ENCSR079VIJ_G401ENCSR401TBQ_Caki2",
  
  #"HepG2",                             # !!! RUNNING !!!
  
  "ENCSR444WCZ_A549",
  "NCI-H460",
  "ENCSR444WCZ_A549NCI-H460",
  
  "K562",
  
  "GSM2334834_U266_HindIII",
  "GSM2334832_RPMI-8226_HindIII",
  "GSM2334834_U266_HindIIIGSM2334832_RPMI-8226_HindIII",
  
  "ENCSR834DXR_SK-N-MC",
  "ENCSR105KFX_SK-N-DZ",
  # NB: no consensus because of metastatic sites
  
  "Panc1_rep12",
  
  "ENCSR346DCU_LNCaP",
  "GSE73782_PC3",              # !!! RUNNING !!!
  "ENCSR346DCU_LNCaPGSE73782_PC3",# !!! RUNNING !!!
  "GSE73782_PC3_ICE",              # !!! RUNNING !!!
  "ENCSR346DCU_LNCaPGSE73782_PC3_ICE",# !!! RUNNING !!!
  
  "ENCSR312KHQ_SK-MEL-5",
  "ENCSR862OGI_RPMI-7951",
  "ENCSR312KHQ_SK-MEL-5ENCSR862OGI_RPMI-7951",
  
  "GSE105194_spinal_cord",
  # "GSE105194_cerebellum"   # !!! RUNNING"
  
  "pipelineConsensus"
)
stopifnot(cl_to_cmp %in% cl_names)
stopifnot(dir.exists(paste0(cl_to_cmp, folderSuffix)))

all_cmps <- combn(cl_to_cmp, m = 2)

stopifnot(nrow(all_cmps) == 2)
stopifnot(ncol(all_cmps) > 0)

i=1
if(buildTable){
all_MoC_dt <- foreach(i = seq_len(ncol(all_cmps)), .combine="rbind") %dopar% {
  
  ds1 <- all_cmps[1,i]
  ds2 <- all_cmps[2,i]
  cat(paste0("*** START: ", ds1, " vs. ", ds2, "\n"))
  
  folder1 <- paste0(ds1, folderSuffix)
  stopifnot(dir.exists(folder1))
  all_files1 <- list.files(file.path(folder1, TADfolder), pattern=TADfilePattern, recursive=FALSE, full.names = TRUE)
  stopifnot(length(all_files1) > 0)  
  stopifnot(length(all_files1) <= 23)
  all_chromos1 <- sapply(all_files1, function(file1) {
    chromo <- gsub(".+_(chr.+?)_.+", "\\1", basename(file1))
    stopifnot(length(chromo) == 1)
    stopifnot(chromo %in% all_chromos)
    chromo
  })
  
  folder2 <- paste0(ds2, folderSuffix)
  stopifnot(dir.exists(folder2))
  all_files2 <- list.files(file.path(folder2, TADfolder), pattern=TADfilePattern, recursive=FALSE, full.names = TRUE)
  stopifnot(length(all_files2) > 0)  
  stopifnot(length(all_files2) <= 23)
  all_chromos2 <- sapply(all_files2, function(file2) {
    chromo <- gsub(".+_(chr.+?)_.+", "\\1", basename(file2))
    stopifnot(length(chromo) == 1)
    stopifnot(chromo %in% all_chromos)
    chromo
  })
  
  intersectChromos <- intersect(all_chromos1, all_chromos2)
  stopifnot(length(intersectChromos) > 0)
  stopifnot(intersectChromos %in% all_chromos)
  
  cat(paste0("... ", ds1, " vs. ", ds2, " - # intersectChromos:", length(intersectChromos), "\n"))
  
  chromo = "chr1"
  chr_MoC_dt <- foreach(chromo = intersectChromos, .combine='rbind') %do% {
    
    cat(paste0("... ", ds1, " vs. ", ds2, " - ", chromo, "\n"))
    
    file1 <- all_files1[grepl(paste0("_", chromo, "_"), basename(all_files1)) & grepl(paste0("^", ds1, "_"), all_files1)]
    stopifnot(length(file1) == 1)
    
    file2 <- all_files2[grepl(paste0("_", chromo, "_"), basename(all_files2)) & grepl(paste0("^", ds2, "_"), all_files2)]
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
    
    # because depending on the processing it could be 1 bp missing at the last domain
    last_before1 <- dt1$end[nrow(dt1)]
    last_after1 <- ceiling(last_before1/binSize)*binSize
    stopifnot(last_after1 >= last_before1)
    dt1$end[nrow(dt1)] <- last_after1
    stopifnot(dt1$end %% binSize == 0)
    stopifnot( (dt1$start-1) %% binSize == 0)
    if(last_before1 != last_after1) {
      txt <- paste0("... change last end from\t", last_before1, "\tto\t", last_after1, "\n" )
      printAndLog(txt, logFile)
    }
    
    last_before2 <- dt2$end[nrow(dt2)]
    last_after2 <- ceiling(last_before2/binSize)*binSize
    stopifnot(last_after2 >= last_before2)
    dt2$end[nrow(dt2)] <- last_after2
    stopifnot(dt2$end %% binSize == 0)
    stopifnot( (dt2$start-1) %% binSize == 0)
    if(last_before2 != last_after2){
      txt <- paste0("... change last end from\t", last_before2, "\tto\t", last_after2, "\n" )
      printAndLog(txt, logFile)
    }
    
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
} else {
  outFile <- file.path(outFold, "all_MoC_dt.Rdata")
  all_MoC_dt <- eval(parse(text = load(outFile)))
}

# 
# all_MoC_dt <- eval(parse(text = load(outFile)))
# stopifnot(nrow(all_MoC_dt) == (length(intersectChromos) * ncol(all_cmps)))
# load("CMP_DATASETS_MOC/all_MoC_dt.Rdata")
# stop("-- ok -- \n")

# take intersect chromos here
tmpDT <- all_MoC_dt
tmpDT$cmp <- paste0(tmpDT$ds1, tmpDT$ds2) 
intersectChromos <- Reduce(intersect, by(tmpDT$cmp, data=tmpDT, FUN=function(x) unique(as.character(x$chromo))))
rm(tmpDT)
txt <- paste0("... # intersectChromos:\t", length(intersectChromos), "\n")
printAndLog(txt, logFile)
txt <- paste0("... intersectChromos:\t", paste0(intersectChromos, collapse=","), "\n")
printAndLog(txt, logFile)

all_MoC_dt$ds1_label <- sapply(as.character(all_MoC_dt$ds1), function(x) {   # !!! NEED THE AS.CHARACTER HERE !!!
  dslab <- as.character(cl_labs[x])
  stopifnot(length(dslab) == 1)
  if(is.na(dslab)) {
    paste0(names(cl_names[cl_names == x]))
  } else {
    paste0(names(cl_names[cl_names == x]), "\n(", dslab, ")")
  }
})
all_MoC_dt$ds2_label <- sapply(as.character(all_MoC_dt$ds2), function(x) {   # !!! NEED THE AS.CHARACTER HERE !!!
  dslab <- as.character(cl_labs[x])
  stopifnot(length(dslab) == 1)
  if(is.na(dslab)) {
    paste0(names(cl_names[cl_names == x]))
  } else {
    paste0(names(cl_names[cl_names == x]), "\n(", dslab, ")")
  }
})

stopifnot(nrow(all_MoC_dt) == (length(intersectChromos) * ncol(all_cmps)))
all_MoC_dt_save <- all_MoC_dt

#******************************************************************************************************************************************** DRAW SYMMETRIC MATRIX

mean_MoC_dt <- aggregate(as.formula(paste0("MoC ~ ds1", xlabType, " + ds2", xlabType)), FUN=mean, data = all_MoC_dt)
stopifnot(!is.na(mean_MoC_dt))

mean_MoC_dt1 <- mean_MoC_dt[,c(paste0("ds1", xlabType), paste0("ds2", xlabType), "MoC")]
mean_MoC_dt2 <- mean_MoC_dt[,c(paste0("ds2", xlabType), paste0("ds1", xlabType), "MoC")]

# self_MoC_dt3 <- data.frame(ds1=all_ds, ds2=all_ds, MoC = 1)
self_MoC_dt3 <- data.frame(ds1=unique(c(all_MoC_dt[, paste0("ds1", xlabType)], all_MoC_dt[, paste0("ds2", xlabType)])),
                           ds2=unique(c(all_MoC_dt[, paste0("ds1", xlabType)], all_MoC_dt[, paste0("ds2", xlabType)])),
                           MoC = 1)

colnames(mean_MoC_dt2) <- colnames(mean_MoC_dt1) <- colnames(self_MoC_dt3) <- c(paste0("ds1", xlabType), paste0("ds2", xlabType), "MoC")
mocDT <- rbind(rbind( mean_MoC_dt1, mean_MoC_dt2), self_MoC_dt3)
stopifnot(!duplicated(mocDT))

mocDT <- mocDT[order(mocDT[, paste0("ds1", xlabType)], mocDT[, paste0("ds2", xlabType)]),]
corMat <- reshape(mocDT, idvar=paste0("ds1", xlabType), timevar=paste0("ds2", xlabType), direction="wide")
rownames(corMat) <- corMat[, paste0("ds1", xlabType)]
colnames(corMat) <- sub("MoC.", "", colnames(corMat))
corMat[, paste0("ds1", xlabType)] <- NULL
corMat <- as.matrix(corMat)
stopifnot(rownames(corMat) == colnames(corMat))
stopifnot(isSymmetric(as.matrix(corMat)))
stopifnot(!is.na(corMat))

tit <- paste0("MoC between tissues (with consensus)\n")

# with pdf output a different unicode character -> save as svg
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
                                                   comparisonName = "tissues",
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

tmpDS <- unique(c(all_MoC_dt$ds1_label, all_MoC_dt$ds2_label))
consensusTissues <- tmpDS[grepl("Consensus", tmpDS)] 

if(length(consensusTissues) > 0) consTissue=consensusTissues[1]

for(consTissue in consensusTissues) {
  
  tissue <- gsub("Consensus", "", consTissue)
  
  if(tissue == "pipeline") {
    # consensus_dt <- all_MoC_dt[ (grepl(paste0(tissue, "Consensus"), all_MoC_dt$ds1_label) | 
    #                                grepl(paste0(tissue, "Consensus"), all_MoC_dt$ds2_label) )
    #                               ,]
    consensus_dt <- all_MoC_dt[ (grepl(paste0("^", consTissue, "$"), all_MoC_dt$ds1_label) | 
                                   grepl(paste0("^", consTissue, "$"), all_MoC_dt$ds2_label) )
                                ,]
    curr_tit <- paste0("MoC with ", tissue, " consensus")
    
  } else {
    # consensus_dt <- all_MoC_dt[ (grepl(paste0(tissue, "Consensus"), all_MoC_dt$ds1_label) | 
    #                                grepl(paste0(tissue, "Consensus"), all_MoC_dt$ds2_label) ) &
    #                                 (grepl(tolower(tissue), tolower(all_MoC_dt$ds1_label)) & 
    #                                    grepl(tolower(tissue), tolower(all_MoC_dt$ds2_label)) )
    #                               ,]
    consensus_dt <- all_MoC_dt[ (grepl(paste0("^", consTissue, "$"), all_MoC_dt$ds1_label) | 
                                   grepl(paste0("^", consTissue, "$"), all_MoC_dt$ds2_label) ) &
                                  (grepl(tolower(tissue), tolower(all_MoC_dt$ds1_label)) & 
                                     grepl(tolower(tissue), tolower(all_MoC_dt$ds2_label)) )
                                ,]
    curr_tit <- paste0("MoC between ", tissue, " cell lines and ", tissue, " consensus")
  }
  
  stopifnot(nrow(consensus_dt) > 0)
  
  # put the consensusDS in newDS1 column
  # consensus_dt$newDS1 <- ifelse(grepl(paste0(tissue, "Consensus"), consensus_dt$ds1_label), 
  #                               consensus_dt[, paste0("ds1", xlabType)], consensus_dt[, paste0("ds2", xlabType)])
  # consensus_dt$newDS2 <- ifelse(grepl(paste0(tissue, "Consensus"), consensus_dt$ds2_label), 
  #                               consensus_dt[, paste0("ds1", xlabType)], consensus_dt[, paste0("ds2", xlabType)])
  consensus_dt$newDS1 <- ifelse(grepl(paste0("^", consTissue, "$"), consensus_dt$ds1_label), 
                                consensus_dt[, paste0("ds1", xlabType)], consensus_dt[, paste0("ds2", xlabType)])
  consensus_dt$newDS2 <- ifelse(grepl(paste0("^", consTissue, "$"), consensus_dt$ds2_label), 
                                consensus_dt[, paste0("ds1", xlabType)], consensus_dt[, paste0("ds2", xlabType)])
  if(xlabType == "") {
    consensus_dt$newDS1_label <- sapply(as.character(consensus_dt$newDS1), function(x) {   # !!! NEED THE AS.CHARACTER HERE !!!
      dslab <- as.character(cl_labs[x])
      stopifnot(length(dslab) == 1)
      if(is.na(dslab)) {
        paste0(names(cl_names[cl_names == x]))
      } else {
        paste0(names(cl_names[cl_names == x]), "\n(", dslab, ")")
      }
    })
    consensus_dt$newDS2_label <- sapply(as.character(consensus_dt$newDS2), function(x) {   # !!! NEED THE AS.CHARACTER HERE !!!
      dslab <- as.character(cl_labs[x])
      stopifnot(length(dslab) == 1)
      if(is.na(dslab)) {
        paste0(names(cl_names[cl_names == x]))
      } else {
        paste0(names(cl_names[cl_names == x]), "\n(", dslab, ")")
      }
    })
    
  } else if(xlabType == "_label") {
    consensus_dt$newDS1_label <- consensus_dt$newDS1
    consensus_dt$newDS2_label <- consensus_dt$newDS2
  }
  # stopifnot(grepl(paste0(tissue, "Consensus"), consensus_dt$newDS1_label))
  # stopifnot(!grepl(paste0(tissue, "Consensus"), consensus_dt$newDS2_label))
  stopifnot(grepl(paste0("^", consTissue, "$"), consensus_dt$newDS1_label))
  stopifnot(!grepl(paste0("^", consTissue, "$"), consensus_dt$newDS2_label))
  

  consensus_dt$comp <- paste0(consensus_dt$newDS1, "_", consensus_dt$newDS2)
  stopifnot(sapply(seq_len(nrow(consensus_dt)),function(i) grepl(strsplit(consensus_dt[, paste0("ds1", xlabType)][i], "\n")[[1]][1],
                                                                 #consensus_dt[, paste0("ds1", xlabType)][i],
                                                                 consensus_dt$comp[i])))
  stopifnot(sapply(seq_len(nrow(consensus_dt)),function(i) grepl(
    strsplit(consensus_dt[, paste0("ds2", xlabType)][i], "\n")[[1]][1],
    consensus_dt[, paste0("ds2", xlabType)][i], 
    consensus_dt$comp[i])))
  
  mean_consensus_dt <- aggregate(as.formula(paste0("MoC ~ newDS1 + newDS2")), FUN=mean, data = consensus_dt)
  mean_consensus_dt <- mean_consensus_dt[order(mean_consensus_dt[, "MoC"], decreasing = TRUE),]
  consensus_dt$newDS2 <- factor(as.character(consensus_dt$newDS2), levels = mean_consensus_dt$newDS2)
  consensus_dt$chromo <- factor(as.character(consensus_dt$chromo), levels = all_chromos)
  
  stopifnot(!is.na(consensus_dt))
  
  p_common <- ggplot(consensus_dt, aes(x = newDS2, y = MoC)) + 
    geom_boxplot(outlier.shape=NA) +
    scale_x_discrete(name="")+
    scale_y_continuous(name=paste0("MoC with consensus"),
                       breaks = scales::pretty_breaks(n = 10))+ 
    labs(colour  = "") +
    ggtitle(label = paste0(curr_tit))+
    theme( # Increase size of axis lines
      # top, right, bottom and left
      # plot.margin = unit(c(1, 1, 4.5, 1), "lines"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size=10),
      panel.grid = element_blank(),
      axis.text.x = element_text( hjust=1,vjust = 0.5, size=12, angle = 90),
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
  
  # outFile <- file.path(outFold, paste0("MoC_", tissue, "_consensus_boxplot_chromoDots.", plotType))
  # ggsave(plot=p_dot, file = outFile, width = widthBoxplot, height = heightBoxplot)
  # cat(paste0("... written: ", outFile, "\n"))
  # foo <- try(dev.off())
  
  outFile <- file.path(outFold, paste0("MoC_matching_", tissue, "_consensus_boxplot_chromoLabs.", plotType))
  ggsave(plot=p_txt, file = outFile, width = widthBoxplot, height = heightBoxplot)
  cat(paste0("... written: ", outFile, "\n"))
  foo <- try(dev.off())
}

######################################################################################
cat(paste0("... written: ", logFile, "\n"))
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))




