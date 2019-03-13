startTime <- Sys.time()

library(foreach)
library(doMC)
library(ggplot2)
library(GenomicRanges)
library(VennDiagram)

options(scipen=100)

cat("> START: cmp_datasets_Venn.R\n")
# Rscript cmp_datasets_Venn.R

buildTable <- TRUE

SSHFS <- FALSE
setDir <- ifelse(SSHFS, "/media/electron", "")
setDir <- ifelse(SSHFS, "~/media/electron", "")

source(file.path("utils_fct.R"))
# source("../Dixon2018_integrative_data/MoC_heatmap_fct.R")
source(file.path("datasets_settings.R"))

# if(SSHFS) setwd("/media/electron/mnt/etemp/marie/Cancer_HiC_data_TAD_DA")
if(SSHFS) setwd("~/media/electron/mnt/etemp/marie/Cancer_HiC_data_TAD_DA")

registerDoMC(ifelse(SSHFS, 2, 40))

outFold <- file.path("CMP_DATASETS_VENN")
system(paste0("mkdir -p ", outFold))

logFile <- file.path(outFold, "cmp_datasets_venn_logFile.txt")
if(!SSHFS) system(paste0("rm -f ", logFile))
if(SSHFS) logFile <- ""

plotType <- "svg"
myHeightV <- ifelse(plotType=="png", 300, 7)
myWidthV <- ifelse(plotType=="png", 300, 7)

dsCol1 <- "dodgerblue4"
dsCol2 <- "darkorange2" 

binSize <- 40000
binSizeKb <- binSize/1000
tolRad <- 2*binSize

txt <- paste0("!! Hard-coded buildTable:\t", as.character(buildTable), "\n")
printAndLog(txt, logFile)
txt <- paste0("!! hard-coded bin size:\t", binSize, "\n")
printAndLog(txt, logFile)
txt <- paste0("!! hard-coded tolerance radius:\t", tolRad, "\n")
printAndLog(txt, logFile)

folderSuffix <- paste0("_", binSizeKb, "kb")
TADfilePattern <- "_final_domains.txt$"
TADfolder <- "FINAL_DOMAINS"

all_chromos <- c(paste0("chr", 1:23), "chrX")  # also used for ordering on the plots

xlabType <- "_label"
stopifnot(xlabType %in% c("", "_label"))


### SELECT WHICH CELL LINES TO INCLUDE IN THE COMPARISONS
# 30 incl. pipCons
cl_to_cmp <- list(
  breast = c("MCF-7",
  "ENCSR549MGQ_T47D",
  "MCF-7ENCSR549MGQ_T47D"),
  
  kidney = c(
  "ENCSR079VIJ_G401",
  "ENCSR401TBQ_Caki2",
  "ENCSR079VIJ_G401ENCSR401TBQ_Caki2"),
  
  lung = c(
  "ENCSR444WCZ_A549",
  "NCI-H460",
  "ENCSR444WCZ_A549NCI-H460"),
  
  myeloma=c(
  "GSM2334834_U266_HindIII",
  "GSM2334832_RPMI-8226_HindIII",
  "GSM2334834_U266_HindIIIGSM2334832_RPMI-8226_HindIII"),
  
  prostate = c(  "ENCSR346DCU_LNCaP",
                 "GSE73782_PC3",              
                 "ENCSR346DCU_LNCaPGSE73782_PC3"),
  prostateICE=c(
  "ENCSR346DCU_LNCaP",
  "GSE73782_PC3_ICE",             
  "ENCSR346DCU_LNCaPGSE73782_PC3_ICE"),
  
  skin=c(
  "ENCSR312KHQ_SK-MEL-5",
  "ENCSR862OGI_RPMI-7951",
  "ENCSR312KHQ_SK-MEL-5ENCSR862OGI_RPMI-7951"),
  
  glioma = c(
  "GSE105194_spinal_cord",
  "GSE105194_cerebellum",  
  "GSE105194_spinal_cordGSE105194_cerebellum")
)

all_ds <- unlist(cl_to_cmp)
cat(all_ds[!all_ds %in% cl_names], "\n")
stopifnot(unlist(lapply(cl_to_cmp, function(x) x[[3]] == paste0(x[[1]], x[[2]]))))
stopifnot(all_ds %in% cl_names)
stopifnot(dir.exists(paste0(all_ds, folderSuffix)))

i=1

allCl_allChr_GRs <- foreach(i = seq_along(cl_to_cmp)) %dopar% {
  
  curr_ds <- cl_to_cmp[[i]]
  curr_folders <- paste0(curr_ds, folderSuffix)
  stopifnot(dir.exists(curr_folders))
  curr_files <- lapply(curr_folders, function(x) {
    list.files(file.path(x, TADfolder), pattern=TADfilePattern, recursive=FALSE, full.names = TRUE)
  })
  names(curr_files) <- curr_folders
  stopifnot(length(unique(unlist(lapply(curr_files, length))))==1)
  
  cl_chromos <- lapply(curr_files, function(all_f) {
    chromo <- gsub(".+_(chr.+?)_.+", "\\1", basename(all_f))
    # stopifnot(length(chromo) == 1)
    stopifnot(chromo %in% all_chromos)
    chromo
  })
  stopifnot(length(unique(cl_chromos)) == 1)
  cl_chromos <- unique(cl_chromos)
  
  cl=curr_files[[1]]

  oneCl_allChr_GRs <- lapply(seq_along(curr_files), function(i_cl){
    
    cl <- curr_files[[i_cl]]
    cl_name <- names(curr_files)[i_cl]
    
    oneCl_allChr_DT <- foreach(clFile = cl, .combine='rbind') %do% {
    dt <- read.delim(clFile, stringsAsFactors = FALSE, header=F, col.names = c("chromo", "start", "end"))
    if(is.character(dt[1,2]) & is.character(dt[2,2])) {
      dt <- read.delim(clFile, stringsAsFactors = FALSE, header=T, col.names = c("chromo", "start", "end"))
    }
    if(!is.numeric(dt[, 2])) cat(clFile, "\n")
    stopifnot(ncol(dt) == 3)
    stopifnot(is.numeric(dt[,2]))
    stopifnot(is.numeric(dt[,3]))
    head(dt, 2)
    # because depending on the processing it could be 1 bp missing at the last domain
    last_before1 <- dt$end[nrow(dt)]
    last_after1 <- ceiling(last_before1/binSize)*binSize
    stopifnot(last_after1 >= last_before1)
    dt$end[nrow(dt)] <- last_after1
    stopifnot(dt$end %% binSize == 0)
    stopifnot( (dt$start-1) %% binSize == 0)
    if(last_before1 != last_after1) {
      txt <- paste0("... change last end from\t", last_before1, "\tto\t", last_after1, "\n" )
      printAndLog(txt, logFile)
    }
    stopifnot(!is.na(dt))
    stopifnot(dt$end > dt$start)
    dt$start <- dt$start - 1
    dt$ID <- paste0(cl_name, "_TAD", 1:nrow(dt))
    dt
    } # end-foreach iterating over all the chromos of a given CELL LINE -> build table
    
    ### !!! EXTEND WITH RADIUS ON LEFT AND RIGHT
    # txt <- paste0("... start and end of the TADs extended by tolRad = ",tolRad,  "\n")
    # printAndLog(txt, logFile)
    # cl_allChr_IR <- IRanges(start = oneCl_allChr_DT$start - tolRad, 
    #                     width = (oneCl_allChr_DT$end - oneCl_allChr_DT$start + 1 + tolRad), 
    #                       names=oneCl_allChr_DT$ID)
    # version 2 -> +/- tolRad done later and separately for boundary overlap vs TAD matching
    cl_allChr_IR <- IRanges(start = oneCl_allChr_DT$start,
                            width = (oneCl_allChr_DT$end - oneCl_allChr_DT$start + 1), 
                            names=oneCl_allChr_DT$ID)
    cl_allChr_GR <- GRanges(ranges = cl_allChr_IR,
                        seqnames=oneCl_allChr_DT$chromo)
    # stop("-ok")
    # cat(clFile, "\n")
    cl_allChr_GR
  }) # end lapply construct GR for 1 cell line and all chromos
  # names(oneCl_allChr_GRs) <- unlist(tmp_chromos)
  
  names(oneCl_allChr_GRs) <- names(curr_files)
  oneCl_allChr_GRs
} # end foreach construct GR for all cell lines and all chromos

names(allCl_allChr_GRs) <- names(cl_to_cmp)

outFile <- file.path(outFold, "allCl_allChr_GRs.Rdata")
save(allCl_allChr_GRs, file = outFile)
cat(paste0("... written: ", outFile, "\n"))

# stop("-- ok --\n")
# load("CMP_DATASETS_VENN/allCl_allChr_GRs.Rdata")

i=1
for(i in seq_along(allCl_allChr_GRs)) {
    
  tissue_name <- names(allCl_allChr_GRs)[i]
  tissue_data <- allCl_allChr_GRs[[i]]
  tissue_cls_cmb <- combn(x = names(tissue_data), m=2)
  #query %over% subject
  # mergeByOverlaps computes the overlap between query and subject according to the arguments in .... 
  # It then extracts the corresponding hits from each object and returns a DataFrame containing one column for the query and one for the subject, 
  # as well as any mcols that were present on either object. The query and subject columns are named by quoting and deparsing the corresponding argument.
  i_col=1
  for(i_col in seq_len(ncol(tissue_cls_cmb))) {
    
    ds1 <- tissue_cls_cmb[1,i_col]
    ds2 <- tissue_cls_cmb[2,i_col]
    GR1 <- tissue_data[[paste0(ds1)]]
    GR2 <- tissue_data[[paste0(ds2)]]
    GR1_2 <- c(GR1, GR2)
    
    ############################################################## STRICT TAD MATCHING
    
    tad_DT <- data.frame(
      start = c(start(GR1), start(GR2)),
      end = c(end(GR1), end(GR2)),
      chromo = c(as.vector(seqnames(GR1)), as.vector(seqnames(GR2))),
      dataset = c(rep(ds1, length(GR1)), rep(ds2, length(GR2))),
      id = c(    paste0( rep(ds1, length(GR1)), as.vector(seqnames(GR1)), "_TAD", 1:length(GR1)),
                 paste0( rep(ds2, length(GR2)), as.vector(seqnames(GR2)), "_TAD", 1:length(GR2))),
      stringsAsFactors = FALSE
    )
  
    i_row=1
    has_match <- foreach(i_row = seq_len(nrow(tad_DT)), .combine = 'c') %dopar% {
      curr_start <- tad_DT$start[i_row]
      curr_end <- tad_DT$end[i_row]
      curr_chromo <- tad_DT$chromo[i_row]
      curr_dataset <- tad_DT$dataset[i_row]
      
      # strict match for start and end ?
      any(abs(curr_start - tad_DT$start) <= tolRad &
        abs(curr_end - tad_DT$end) <= tolRad & 
        curr_chromo == tad_DT$chromo & curr_dataset != tad_DT$dataset)
    }
    nMatch_ds1 <- sum(as.numeric(has_match & tad_DT$dataset == ds1))
    nMatch_ds1
    nMatch_ds2 <- sum(as.numeric(has_match & tad_DT$dataset == ds2))
    nMatch_ds2
    # stopifnot(nMatch_ds1 == nMatch_ds2)
    intersect_ds1_ds2 <- nMatch_ds1
      
    uniq_ds1 <- sum(as.numeric(!has_match & tad_DT$dataset == ds1))
    uniq_ds2 <- sum(as.numeric(!has_match & tad_DT$dataset == ds2))
    
    txt <- paste0("nMatch_ds1 = ", nMatch_ds1, "\n")
    printAndLog(txt, logFile)
    txt <- paste0("nMatch_ds2 = ", nMatch_ds2, "\n")
    printAndLog(txt, logFile)
    # stopifnot(nMatch_ds1 == nMatch_ds2)
    txt <- paste0("intersect_ds1_ds2 = ", intersect_ds1_ds2, "\n")
    printAndLog(txt, logFile)
    txt <- paste0("uniq_ds1 = ", uniq_ds1, "\n")
    printAndLog(txt, logFile)
    txt <- paste0("uniq_ds2 = ", uniq_ds2, "\n")
    printAndLog(txt, logFile)
    
    outFile <- file.path(outFold, paste0(ds1, "_", ds2, "_pairwiseVenn_tolRad_", tolRad, "_strictTAD.", plotType))
    do.call(plotType, list(outFile, height=myHeightV, width=myWidthV))
    grid.newpage()
    draw.pairwise.venn(area1 = uniq_ds1+intersect_ds1_ds2, 
                       area2 = uniq_ds2 + intersect_ds1_ds2, 
                       cross.area =  intersect_ds1_ds2,
                       category = c(ds1, ds2), 
                       lty = "blank", 
                       alpha = rep(0.5, 2),
                       fill = c(dsCol1, dsCol2))
    grid.text(y = 0.9, paste0(tissue_name, ": strict TAD matching"), 
              just=c("centre", "top"), gp=gpar(fontsize=16,font=2))
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
  
    ############################################################## BOUNDARY MATCHING
    bd_start_DT <- tad_DT
    bd_start_DT$end <- bd_start_DT$start
    colnames(bd_start_DT)[colnames(bd_start_DT) == "start"] <- "bd_left"
    colnames(bd_start_DT)[colnames(bd_start_DT) == "end"] <- "bd_right"
    bd_start_DT$bd_left <- bd_start_DT$bd_left - tolRad
    bd_start_DT$bd_right <- bd_start_DT$bd_right + tolRad
    bd_start_DT$pos <- "start"
    
    bd_end_DT <- tad_DT
    bd_end_DT$start <- bd_end_DT$end
    colnames(bd_end_DT)[colnames(bd_end_DT) == "start"] <- "bd_left"
    colnames(bd_end_DT)[colnames(bd_end_DT) == "end"] <- "bd_right"
    bd_end_DT$bd_left <- bd_end_DT$bd_left - tolRad
    bd_end_DT$bd_right <- bd_end_DT$bd_right + tolRad
    bd_end_DT$pos <- "end"
    
    all_bd_DT <- rbind(bd_start_DT, bd_end_DT)
    
    nTotBD <- nrow(all_bd_DT)
    txt <- paste0("nTotBD = ", nTotBD, "\n")
    printAndLog(txt, logFile)
    
    all_bd_IR <- IRanges(start = all_bd_DT$bd_left,
                            width = (all_bd_DT$bd_right - all_bd_DT$bd_left + 1), 
                            names=all_bd_DT$dataset)
    all_bd_GR <- GRanges(ranges = all_bd_IR,
                            seqnames=all_bd_DT$chromo)
    
    nTotBD_GR <- length(all_bd_GR)
    txt <- paste0("nTotBD_GR = ", nTotBD_GR, "\n")
    printAndLog(txt, logFile)
    stopifnot(nTotBD == nTotBD_GR)
    
    start(all_bd_GR) <- start(all_bd_GR) + 1
    
    ds1BD <- all_bd_GR[names(all_bd_GR) == ds1]  %over% all_bd_GR[names(all_bd_GR) == ds2]  
    ds2BD <- all_bd_GR[names(all_bd_GR) == ds2]  %over% all_bd_GR[names(all_bd_GR) == ds1]
    
    nMatchBD_ds1 <- sum(ds1BD)
    nMatchBD_ds2 <- sum(ds2BD)
    
    noMatchBD_ds1 <- sum(!ds1BD)
    noMatchBD_ds2 <- sum(!ds2BD)
    
    nBD_ds1 <- length(ds1BD)
    nBD_ds2 <- length(ds2BD)
    
    txt <- paste0("nMatchBD_ds1 = ", nMatchBD_ds1, "\n")
    printAndLog(txt, logFile)
    txt <- paste0("nMatchBD_ds2 = ", nMatchBD_ds2, "\n")
    printAndLog(txt, logFile)
    txt <- paste0("noMatchBD_ds1 = ", noMatchBD_ds1, "\n")
    printAndLog(txt, logFile)
    txt <- paste0("noMatchBD_ds2 = ", noMatchBD_ds2, "\n")
    printAndLog(txt, logFile)
    txt <- paste0("nBD_ds1 = ", nBD_ds1, "\n")
    printAndLog(txt, logFile)
    txt <- paste0("nBD_ds2 = ", nBD_ds2, "\n")
    printAndLog(txt, logFile)
    
    dup_bd_GR <- all_bd_GR[all_bd_GR %in% all_bd_GR[duplicated(all_bd_GR)]]
    uniq_bd_GR <- all_bd_GR[!all_bd_GR %in% all_bd_GR[duplicated(all_bd_GR)]]
    nTot_allBD <- length(all_bd_GR) # 20006
    txt <- paste0("nTot_allBD = ", nTot_allBD, "\n")
    printAndLog(txt, logFile)
    
    nTot_allBD_uniq <- length(uniq_bd_GR)  # 11573
    txt <- paste0("nTot_allBD_uniq = ", nTot_allBD_uniq, "\n")
    printAndLog(txt, logFile)
    
    # length(unique(all_bd_GR)) # 15663
    # length(dup_bd_GR) # 8433
    # sum(duplicated(all_bd_GR)) # 4343    
    # # 20006-15663
    # [1] 4343
    
    # reduce first orders the ranges in x from left to right, then merges the overlapping or adjacent ones.
    all_bd_GR_rd <- reduce(all_bd_GR)
    nTotBD_GR_rd <- length(all_bd_GR_rd)
    txt <- paste0("nTotBD_GR_rd = ", nTotBD_GR_rd, "\n")
    printAndLog(txt, logFile)
    
    all_bd_GR_uniq <- unique(all_bd_GR)
    nTotBD_GR_uniq <- length(all_bd_GR_uniq)
    txt <- paste0("nTotBD_GR_uniq = ", nTotBD_GR_rd, "\n")
    printAndLog(txt, logFile)
    
    outFile <- file.path(outFold, paste0(ds1, "_", ds2, "_pairwiseVenn_matchBD.", plotType))
    do.call(plotType, list(outFile, height=myHeightV, width=myWidthV))
    grid.newpage()
    draw.pairwise.venn(area1 = noMatchBD_ds1+intersect_ds1_ds2, 
                       area2 = noMatchBD_ds2 + intersect_ds1_ds2, 
                       cross.area =  intersect_ds1_ds2,
                       category = c(ds1, ds2), 
                       lty = "blank", 
                       alpha = rep(0.5, 2),
                       fill = c(dsCol1, dsCol2))
    grid.text(y = 0.9, paste0(tissue_name, ": boundary matching"), 
              just=c("centre", "top"), gp=gpar(fontsize=16,font=2))
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
  } # end-for iterating over pair of cell line
} # end-for iterating over tissue


######################################################################################
cat(paste0("... written: ", logFile, "\n"))
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))






