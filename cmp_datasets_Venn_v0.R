startTime <- Sys.time()

library(foreach)
library(doMC)
library(ggplot2)
library(GenomicRanges)

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
widthMat <- 26
heightMat <- 14
widthBoxplot <- 10
heightBoxplot <- 7

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
    cat(clFile, "\n")
    cl_allChr_GR
  }) # end lapply construct GR for 1 cell line and all chromos
  # names(oneCl_allChr_GRs) <- unlist(tmp_chromos)
  
  names(oneCl_allChr_GRs) <- names(curr_files)
  oneCl_allChr_GRs
} # end foreach construct GR for all cell lines and all chromos

names(allCl_allChr_GRs) <- cl_to_cmp

outFile <- file.path(outFold, "allCl_allChr_GRs.Rdata")
save(allCl_allChr_GRs, file = outFile)
cat(paste0("... written: ", outFile, "\n"))

stop("-- ok --\n")

load("CMP_DATASETS_VENN/allCl_allChr_GRs.Rdata")


i=1

tissue_data <- allCl_allChr_GRs[[1]]

#query %over% subject

# mergeByOverlaps computes the overlap between query and subject according to the arguments in .... 
# It then extracts the corresponding hits from each object and returns a DataFrame containing one column for the query and one for the subject, 
# as well as any mcols that were present on either object. The query and subject columns are named by quoting and deparsing the corresponding argument.

tissue_cls_cmb <- combn(x = names(tissue_data), m=2)

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
  stopifnot(nMatch_ds1 == nMatch_ds2)
  intersect_ds1_ds2 <- nMatch_ds1
    
  uniq_ds1 <- sum(as.numeric(!has_match & tad_DT$dataset == ds1))
  uniq_ds2 <- sum(as.numeric(!has_match & tad_DT$dataset == ds2))
  
  
  
  txt <- paste0("nMatch_ds1 = ", nMatch_ds1, "\n")
  cat(txt, logFile)
  txt <- paste0("nMatch_ds2 = ", nMatch_ds2, "\n")
  cat(txt, logFile)
  stopifnot(nMatch_ds1 == nMatch_ds2)
  txt <- paste0("intersect_ds1_ds2 = ", intersect_ds1_ds2, "\n")
  cat(txt, logFile)
  txt <- paste0("uniq_ds1 = ", uniq_ds1, "\n")
  cat(txt, logFile)
  txt <- paste0("uniq_ds2 = ", uniq_ds2, "\n")
  cat(txt, logFile)

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
  cat(txt, logFile)
  
  all_bd_IR <- IRanges(start = all_bd_DT$bd_left,
                          width = (all_bd_DT$bd_right - all_bd_DT$bd_left + 1), 
                          names=all_bd_DT$dataset)
  all_bd_GR <- GRanges(ranges = all_bd_IR,
                          seqnames=all_bd_DT$chromo)
  
  nTotBD_GR <- length(all_bd_GR)
  txt <- paste0("nTotBD_GR = ", nTotBD_GR, "\n")
  cat(txt, logFile)
  stopifnot(nTotBD == nTotBD_GR)
  
  start(all_bd_GR) <- start(all_bd_GR) + 1
  
  ds1BD <- all_bd_GR[names(all_bd_GR) == ds1]  %over% all_bd_GR[names(all_bd_GR) == ds2]  
  ds2BD <- all_bd_GR[names(all_bd_GR) == ds2]  %over% all_bd_GR[names(all_bd_GR) == ds1]
  
  nMatchBD_ds1 <- sum(ds1BD)
  nMatchBD_ds2 <- sum(ds2BD)
  
  noMatchBD_ds1 <- sum(!ds1BD)
  noMatchBD_ds2 <- sum(!ds2BD)
  
  # reduce first orders the ranges in x from left to right, then merges the overlapping or adjacent ones.
  all_bd_GR_rd <- reduce(all_bd_GR)
  nTotBD_GR_rd <- length(all_bd_GR_rd)
  txt <- paste0("nTotBD_GR_rd = ", nTotBD_GR_rd, "\n")
  cat(txt, logFile)
  
  all_bd_GR_uniq <- unique(all_bd_GR)
  nTotBD_GR_uniq <- length(all_bd_GR_uniq)
  txt <- paste0("nTotBD_GR_uniq = ", nTotBD_GR_rd, "\n")
  cat(txt, logFile)
  
  
  
  
  
  # otherwise ...-40000 will match with 40000-
  start(GR1_2) <- start(GR1_2) + 1
  
  # reduce first orders the ranges in x from left to right, then merges the overlapping or adjacent ones.
     # drop.empty.ranges	-> TRUE or FALSE. Should empty ranges be dropped?
     # min.gapwidth   Ranges separated by a gap of at least min.gapwidth positions are not merged.
     # with.revmap	  TRUE or FALSE. Should the mapping from output to input ranges be stored in the returned object? If yes, then it is stored as metadata column revmap of type IntegerList.
     # with.inframe.attrib	    TRUE or FALSE. For internal use
    
  
  GR1_2_rd <- reduce(GR1_2)
  
  ds1_overBD <- GR1 %over% GR1_2_rd
  ds2_overBD <- GR2 %over% GR1_2_rd
  
  uniqBD_ds1 <- sum(!ds1_overBD)
  commonBD_ds1 <- sum(ds1_overBD)
  uniqBD_ds2 <- sum(!ds2_overBD)
  commonBD_ds2 <- sum(ds2_overBD)
  
  
  txt <- paste0("GR1_2_rd = ", length(GR1_2_rd), "\n")
  cat(txt, logFile)
  
  txt <- paste0("uniqBD_ds1 = ", uniqBD_ds1, "\n")
  cat(txt, logFile)
  txt <- paste0("uniqBD_ds2 = ", uniqBD_ds2, "\n")
  cat(txt, logFile)  
  txt <- paste0("commonBD_ds1 = ", commonBD_ds1, "\n")
  cat(txt, logFile)
  txt <- paste0("commonBD_ds2 = ", commonBD_ds2, "\n")
  cat(txt, logFile)
  ## => need to convert from TAD to BD !!!
  
  q1_s2 <- 
  stopifnot(length(q1_s2) == length(GR1))
  q2_s1 <- GR2 %over% GR1
  stopifnot(length(q2_s1) == length(GR2))
  GR1_2 <- unique(c(GR1, GR2))
  GR1_2_matching <- mergeByOverlaps(query=GR1_2, GR1_2)
  
  q1_s1_2 <- GR1 %over% GR1_2
  stopifnot(length(q1_s1_2) == length(GR1))
  
  q2_s1_2 <- GR2 %over% GR1_2
  stopifnot(length(q2_s1_2) == length(GR2))
  
  # query %over% subject
  # `%over%` <- function(query, subject) overlapsAny(query, subject) => finds the ranges in query that overlap any of the ranges in subject
  
  merge_ds1_ds2 <- mergeByOverlaps(query=GR1, subject=GR2)
  merge_ds2_ds1 <- mergeByOverlaps(query=GR2, subject=GR1)
  
  merge_ds2_ds1 <- mergeByOverlaps(query=GR1_2, subject=GR1_2)
  
  uniqDS1 <- sum(!q1_s2)
  uniqDS2 <- sum(!q2_s1)
  commonDS1 <- sum(q1_s2)
  commonDS2 <- sum(q2_s1)
  
  uniqDS1_over1_2 <- sum(!q1_s1_2)
  uniqDS2_over1_2 <- sum(!q2_s1_2)
  commonDS1_over1_2 <- sum(q1_s1_2)
  commonDS2_over1_2 <- sum(q2_s1_2)
  
  
  
  txt <- paste0("n TADs DS1 = ", length(GR1), "\n")
  cat(txt, logFile)
  txt <- paste0("n TADs DS2 = ", length(GR2), "\n")
  cat(txt, logFile)
  txt <- paste0("n TADs DS1_2 = ", length(GR1_2), "\n")
  cat(txt, logFile)
  
  
  
  
  
  
  
  
  
  subjectHits(findOverlaps(bdRange, reduce(bdRange)))
  
  
  
  range1_DT <- data.frame(set="set1", start = set1_BD-tolRad, end = set1_BD+tolRad)
  range2_DT <- data.frame(set="set2", start = set2_BD-tolRad, end = set2_BD+tolRad)
  rangeDT <- rbind(range1_DT, range2_DT)
  # IRanges with all starts and ends of both callers
  bdRange <-   IRanges(rangeDT$start, rangeDT$end)
  # reduce first orders the ranges in x from left to right, then merges the overlapping or adjacent ones.
  # the look for matching between reduced boundary regions of both callers with the boundary regions of each caller
  # => because rangeDT is in the same order as bdRange, for each boundary region of rangeDT, get the matching with the merged union of boundary regions
  rangeDT$group <- subjectHits(findOverlaps(bdRange, reduce(bdRange)))
  # for the DT with all boundary regions -> retrieve to which merged union boundary region they have match 
  # then for each of the merged boundary region (x$group), look if this boundary regions has match with regions in set1, set1+set2, set2
  matchDT <- data.frame(do.call(rbind, by(rangeDT, rangeDT$group, function(x) c(unique(x$group), mean(as.numeric(unique(x$set) == "set1"))))))
  colnames(matchDT) <- c("group", "set1_match")
  
  stopifnot(all(matchDT$set1_match %in% c(0,0.5,1)))
  
  nUniqueSet1 <- sum(matchDT$set1_match == 1)
  nUniqueSet2 <- sum(matchDT$set1_match == 0)
  nShared <-  sum(matchDT$set1_match == 0.5)
  
  if(verbose) cat(paste0("... number of unique boundaries in ", set1_name, ": ", n
  
                         
                         
                         
   q1_s2 <- GR1 %over% GR2
   stopifnot(length(q1_s2) == length(GR1))
   q2_s1 <- GR2 %over% GR1
   stopifnot(length(q2_s1) == length(GR2))
   GR1_2 <- unique(c(GR1, GR2))
   GR1_2_matching <- mergeByOverlaps(query=GR1_2, GR1_2)
              
  q1_s1_2 <- GR1 %over% GR1_2
  stopifnot(length(q1_s1_2) == length(GR1))
  
  q2_s1_2 <- GR2 %over% GR1_2
  stopifnot(length(q2_s1_2) == length(GR2))
  
  # query %over% subject
  # `%over%` <- function(query, subject) overlapsAny(query, subject) => finds the ranges in query that overlap any of the ranges in subject
  
  merge_ds1_ds2 <- mergeByOverlaps(query=GR1, subject=GR2)
  merge_ds2_ds1 <- mergeByOverlaps(query=GR2, subject=GR1)
  
  merge_ds2_ds1 <- mergeByOverlaps(query=GR1_2, subject=GR1_2)
  
  uniqDS1 <- sum(!q1_s2)
  uniqDS2 <- sum(!q2_s1)
  commonDS1 <- sum(q1_s2)
  commonDS2 <- sum(q2_s1)
  
  uniqDS1_over1_2 <- sum(!q1_s1_2)
  uniqDS2_over1_2 <- sum(!q2_s1_2)
  commonDS1_over1_2 <- sum(q1_s1_2)
  commonDS2_over1_2 <- sum(q2_s1_2)
  
  
  
  txt <- paste0("n TADs DS1 = ", length(GR1), "\n")
  cat(txt, logFile)
  txt <- paste0("n TADs DS2 = ", length(GR2), "\n")
  cat(txt, logFile)
  txt <- paste0("n TADs DS1_2 = ", length(GR1_2), "\n")
  cat(txt, logFile)

  
  
  txt <- paste0("n unique TADs DS1 = ", uniqDS1, "\n")
  cat(txt, logFile)
  txt <- paste0("n common TADs DS1 = ", commonDS1, "\n")  
  cat(txt, logFile)
  txt <- paste0("n unique TADs DS1 = ", uniqDS1_over1_2, "\n")
  cat(txt, logFile)
  txt <- paste0("n common TADs DS1 over DS1_2  = ", commonDS1_over1_2, "\n")
  cat(txt, logFile)
  
  txt <- paste0("n unique TADs DS2 = ", uniqDS2, "\n")
  cat(txt, logFile)
  txt <- paste0("n common TADs DS2 = ", commonDS2, "\n")  
  cat(txt, logFile)
  txt <- paste0("n unique TADs DS2 over DS1_2 = ", uniqDS2_over1_2, "\n")
  cat(txt, logFile)
  txt <- paste0("n common TADs DS2 = ", commonDS2_over1_2, "\n")  
  cat(txt, logFile)
  
  
  
  
  # stopifnot(commonDS1 == commonDS2)
  
  library(VennDiagram)
  
  grid.newpage()
  draw.triple.venn(area1 = 22, area2 = 20, area3 = 13, n12 = 11, n23 = 4, n13 = 5, 
                   n123 = 1, category = c("Dog People", "Cat People", "Lizard People"), lty = "blank", 
                   fill = c("skyblue", "pink1", "mediumorchid"))
  
  grid.newpage()
  draw.pairwise.venn(area1 = uniqDS1, area2 = uniqDS2, cross.area =  min(commonDS1, commonDS2), 
                 category = c("uniqDS1", "uniqDS2"), lty = "blank", 
                   fill = c("skyblue", "pink1"))
  
  
  grid.newpage()
  grid.newpage()
  draw.pairwise.venn(area1 = uniqDS1, area2 = uniqDS2, cross.area =  0, 
                     category = c("uniqDS1", "uniqDS2"), lty = "blank", 
                     fill = c("skyblue", "pink1"))
  
  
  draw..venn(area1 = uniqDS1, area2 = uniqDS2, n12 = min(commonDS1, commonDS2), n13=0), 
                   category = c("uniqDS1", "uniqDS2", "min12"), lty = "blank", 
                   fill = c("skyblue", "pink1", "mediumorchid"))
  
  
}


cl_data <- tissue_data[[1]]

cl_chr_data <- cl_data[[1]]

stopifnot(unique(as.character(levels(cl_chr_data@seqnames))

curr_chromo <- "chr1"

















all_cmps1 <- combn(cl_to_cmp, m = 2)
all_cmps2 <- all_cmps1[c(2,1),]
all_cmps <- cbind(all_cmps1, all_cmps2) # because asymmetric matching !

stopifnot(nrow(all_cmps) == 2)
stopifnot(ncol(all_cmps) > 0)

i=1
if(buildTable){
all_match_dt <- foreach(i = seq_len(ncol(all_cmps)), .combine="rbind") %dopar% {
  # percent matching of the domains from ds1 in ds2
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
  chr_match_dt <- foreach(chromo = intersectChromos, .combine='rbind') %do% {
    
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
    

    
    strict_domainMatching <- sapply(seq_len(nrow(dt1)), function(k){
      curr_start <- dt1$start[k] 
      curr_end <- dt1$end[k]
      
      matchStart <- which(abs(dt2$start - curr_start) <= tolRad)
      matchEnd <- which(abs(dt2$end - curr_end) <= tolRad)
      
      if(length(matchStart) == 0 | length(matchEnd) == 0) return(0)
      
      if(length(intersect(matchStart, matchEnd)) == 0){
        return(0)
      } else {
        return(1)
      }
    })
    # strict_domainMatching <- unlist(strict_domainMatching) # not needed if 1 value returned [and this should be the case !]
    stopifnot(length(strict_domainMatching) == nrow(dt1))
    head(strict_domainMatching)
    cat("sum = ", sum(unlist(strict_domainMatching)), "\n")
    strictMatchRatio <- sum(strict_domainMatching)/length(strict_domainMatching)
    stopifnot(strictMatchRatio >= 0 & strictMatchRatio <= 1)
    
    loose_domainMatching <- sapply(seq_len(nrow(dt1)), function(k){
      curr_start <- dt1$start[k]
      curr_end <- dt1$end[k]
      
      matchStart <- which(abs( dt2$start - curr_start) <= tolRad)
      matchEnd <- which(abs(dt2$end - curr_end) <= tolRad)
      
      if(length(matchStart) == 0 | length(matchEnd) == 0){
        return(0)
      } else {
        return(1)
      }
    })
    # loose_domainMatching <- unlist(loose_domainMatching)
    stopifnot(length(loose_domainMatching) == nrow(dt1))
    looseMatchRatio <- sum(loose_domainMatching)/length(loose_domainMatching)
    stopifnot(looseMatchRatio >= 0 & looseMatchRatio <= 1)
    
    bd_dt1 <- data.frame(chromo = rep(dt1$chromo,2),
                         position = c(dt1$start, dt1$end),
                         stringsAsFactors = FALSE)
    bd_dt1 <- unique(bd_dt1)
    
    bdMatching <- sapply(seq_len(nrow(bd_dt1)), function(k){
      curr_bd <- bd_dt1$position[k]
      
      matchStart <- which(abs(dt2$start - curr_bd) <= tolRad)
      matchEnd <- which(abs(dt2$end - curr_bd) <= tolRad)
      
      if(length(matchStart) > 0 | length(matchEnd) > 0){
        return(1)
      } else {
        return(0)
      }
    })
    # bdMatching <- unlist(bdMatching)
    stopifnot(length(bdMatching) == nrow(bd_dt1))
    bdMatchRatio <- sum(bdMatching)/length(bdMatching)
    stopifnot(bdMatchRatio >= 0 & bdMatchRatio <= 1)
    
    data.frame(
      ds1 = ds1,
      ds2 = ds2,
      chromo = chromo,
      strictMatchRatio = strictMatchRatio,
      looseMatchRatio = looseMatchRatio,
      bdMatchRatio = bdMatchRatio,
      stringsAsFactors = FALSE
    )
  }
  chr_match_dt
}

outFile <- file.path(outFold, "all_match_dt.Rdata")
save(all_match_dt, file = outFile)
cat(paste0("... written: ", outFile, "\n"))
} else {
  outFile <- file.path(outFold, "all_match_dt.Rdata")
  all_match_dt <- eval(parse(text = load(outFile)))
}
# 
# all_match_dt <- eval(parse(text = load(outFile)))
# stopifnot(nrow(all_match_dt) == (length(intersectChromos) * ncol(all_cmps)))
# load("CMP_DATASETS_MATCHING/all_match_dt.Rdata")
# stop("-- ok -- \n")

# take intersect chromos here
tmpDT <- all_match_dt
tmpDT$cmp <- paste0(tmpDT$ds1, tmpDT$ds2) 
intersectChromos <- Reduce(intersect, by(tmpDT$cmp, data=tmpDT, FUN=function(x) unique(as.character(x$chromo))))
rm(tmpDT)
txt <- paste0("... # intersectChromos:\t", length(intersectChromos), "\n")
printAndLog(txt, logFile)
txt <- paste0("... intersectChromos:\t", paste0(intersectChromos, collapse=","), "\n")
printAndLog(txt, logFile)

all_match_dt$ds1_label <- sapply(as.character(all_match_dt$ds1), function(x) {   # !!! NEED THE AS.CHARACTER HERE !!!
  dslab <- as.character(cl_labs[x])
  stopifnot(length(dslab) == 1)
  if(is.na(dslab)) {
    paste0(names(cl_names[cl_names == x]))
  } else {
    paste0(names(cl_names[cl_names == x]), "\n(", dslab, ")")
  }
})
all_match_dt$ds2_label <- sapply(as.character(all_match_dt$ds2), function(x) {   # !!! NEED THE AS.CHARACTER HERE !!!
  dslab <- as.character(cl_labs[x])
  stopifnot(length(dslab) == 1)
  if(is.na(dslab)) {
    paste0(names(cl_names[cl_names == x]))
  } else {
    paste0(names(cl_names[cl_names == x]), "\n(", dslab, ")")
  }
})
#******************************************************************************************************************************************** DRAW SYMMETRIC MATRIX

var_to_plot <- colnames(all_match_dt)[!colnames(all_match_dt) %in% c("ds1", "ds2", "chromo", "ds1_label", "ds2_label")]

curr_var <- "strictMatchRatio"

plot_tit <- c(
strictMatchRatio = "Strict matching ratio",
looseMatchRatio = "Loose matching ratio",
bdMatchRatio = "Boundary matching ratio"
)

stopifnot(var_to_plot %in% names(plot_tit))

curr_var=var_to_plot[1]
for(curr_var in var_to_plot) {

  mytit <- plot_tit[curr_var]
  
  mean_match_dt <- aggregate(as.formula(paste0(curr_var, " ~ ds1", xlabType, " + ds2", xlabType)), FUN=mean, data = all_match_dt)
  stopifnot(!is.na(mean_match_dt))
  
  self_match_dt <- data.frame(ds1=unique(c(all_match_dt[, paste0("ds1", xlabType)], all_match_dt[, paste0("ds2", xlabType)])), ds2=unique(c(all_match_dt[, paste0("ds1", xlabType)], all_match_dt[, paste0("ds2", xlabType)])), tmp = 1)
  colnames(self_match_dt)[colnames(self_match_dt) == "tmp"] <- curr_var
  colnames(self_match_dt)[colnames(self_match_dt) == "ds1"] <- paste0("ds1", xlabType)
  colnames(self_match_dt)[colnames(self_match_dt) == "ds2"] <- paste0("ds2", xlabType)
  
  ratioDT <- rbind(mean_match_dt, self_match_dt)
  stopifnot(!duplicated(ratioDT))
  
  ratioDT <- ratioDT[order(ratioDT[, paste0("ds1", xlabType)], ratioDT[, paste0("ds2", xlabType)]),]
  corMat <- reshape(ratioDT, idvar=paste0("ds1", xlabType), timevar=paste0("ds2", xlabType), direction="wide")
  rownames(corMat) <- corMat[, paste0("ds1", xlabType)]
  colnames(corMat) <- sub(paste0(curr_var, "."), "", colnames(corMat))
  corMat[, paste0("ds1", xlabType)] <- NULL
  corMat <- as.matrix(corMat)
  stopifnot(rownames(corMat) == colnames(corMat))
  # stopifnot(isSymmetric(as.matrix(corMat)))
  stopifnot(!is.na(corMat))
  corMat[1:3,1:3]
  
  if(any(grepl("Consensus", rownames(corMat)))) {
    tit <- paste0(curr_var, " between tissues (with consensus)\n")
    outFile <- file.path(outFold, paste0(curr_var, "_tissues_with_consensus_match_heatmap.", plotType))
  } else {
    tit <- paste0(curr_var, " between tissues")
    outFile <- file.path(outFold, paste0(curr_var, "_tissues_match_heatmap.", plotType))
  }
  
  outfile_dendro <- file.path(outFold, paste0(curr_var, "row.dend_check.png"))
  # with pdf output a different unicode character -> save as svg
  # outFile <- paste0(outFold, "/", "figure4_match_square_heatmap_with_dendro_", curr_norm, "_", res, "kb_all_callers.", "svg")
  gplot_dendro <- plot_ggheatmap_with_left_rowdendro(x=as.matrix(corMat),
                                                     ranked_branches =T,
                                                     plotMap = "square", 
                                                     low_limit_col = 0,
                                                     high_limit_col = 1,
                                                     fill_legName = paste0(curr_var), 
                                                     dendroLabSize = 4,
                                                     addClusterDot = F,
                                                     annotateMat = TRUE,
                                                     annotateMean = TRUE,
                                                     comparisonName = "caller",
                                                     legCategoryCols = NULL,
                                                     lab_color_vect = NULL)
  
  
  ggsave(plot=gplot_dendro, file = outFile, width = widthMat, height = heightMat)
  cat(paste0("... written: ", outFile, "\n"))
  foo <- try(dev.off())
  
  # outFile <- file.path(outFold, "cmp_match_matrix.pdf")
  # ggsave(plot=gplot_dendro, file = outFile, width = 26, height = 14)
  # cat(paste0("... written: ", outFile, "\n"))
  
  #******************************************************************************************************************************************** BOXPLOT FOR THE CONSENSUS
  tmpDS <- unique(c(all_match_dt$ds1_label, all_match_dt$ds2_label))
  stopifnot(nrow(tmpDS) > 0)
  consensusTissues <- tmpDS[grepl("Consensus", tmpDS)] 
  
  if(length(consensusTissues) > 0) consTissue=consensusTissues[1]
  
  consTissue = "prostateConsensusICE" 
  consTissue = "prostateConsensus" 
  
  for(consTissue in consensusTissues) {
      
    # should also match prostateConsensusICE
    tissue <- gsub("Consensus.*$", "", consTissue)
    
    if(tissue == "pipeline") {
      # consensus_dt <- all_match_dt[ (grepl(paste0(tissue, "Consensus"), all_match_dt$ds1_label) |
      #                                  grepl(paste0(tissue, "Consensus"), all_match_dt$ds2_label) ) 
      #                               ,]
      consensus_dt <- all_match_dt[ (grepl(paste0("^", consTissue, "$"), all_match_dt$ds1_label) |
                                       grepl(paste0("^", consTissue, "$"), all_match_dt$ds2_label) )
                                    ,]
      curr_tit <- paste0(mytit, " with ", tissue, " consensus")
      
    }else {
      # consensus_dt <- all_match_dt[ (grepl(paste0(tissue, "Consensus"), all_match_dt$ds1_label) | 
      #                                  grepl(paste0(tissue, "Consensus"), all_match_dt$ds2_label) ) &
      #                                 (grepl(tolower(tissue), tolower(all_match_dt$ds1_label)) & 
      #                                    grepl(tolower(tissue), tolower(all_match_dt$ds2_label)) )
      #                               ,]
      consensus_dt <- all_match_dt[ (grepl(paste0("^", consTissue, "$"), all_match_dt$ds1_label) | 
                                       grepl(paste0("^", consTissue, "$"), all_match_dt$ds2_label) ) &
                                      (grepl(tolower(tissue), tolower(all_match_dt$ds1_label)) & 
                                         grepl(tolower(tissue), tolower(all_match_dt$ds2_label)) )
                                    ,]
      curr_tit <- paste0(mytit, " between ", tissue, " cell lines and ", tissue, " consensus")
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
    # stopifnot(!grepl(paste0(tissue, "Consensus$"), consensus_dt$newDS2_label))    
    
    cat("*** ", consTissue, "\n")
    
    cat(unique(consensus_dt$newDS1_label), "\n")
    
    cat(unique(consensus_dt$newDS2_label), "\n")
    
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
    
    mean_consensus_dt <- aggregate(as.formula(paste0(curr_var, " ~ newDS1 + newDS2")), FUN=mean, data = consensus_dt)
    mean_consensus_dt <- mean_consensus_dt[order(mean_consensus_dt[, curr_var], decreasing = TRUE),]
    consensus_dt$newDS2 <- factor(as.character(consensus_dt$newDS2), levels = mean_consensus_dt$newDS2)
    consensus_dt$chromo <- factor(as.character(consensus_dt$chromo), levels = all_chromos)
    
    stopifnot(!is.na(consensus_dt))
    
    p_common <- ggplot(consensus_dt, aes_string(x = "newDS2", y = curr_var)) + 
      geom_boxplot(outlier.shape=NA) +
      scale_x_discrete(name="") +
      scale_y_continuous(name=paste0(mytit, " with consensus"),
                         breaks = scales::pretty_breaks(n = 10))+ 
      labs(colour  = "") +
      ggtitle(label = paste0(curr_tit))+
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
      ) 
    
    if(SSHFS) p_common
    p_dot <- p_common + geom_jitter(aes(colour = chromo)) 
    if(SSHFS) p_dot
    p_txt <- p_common + geom_text(aes(label=chromo, colour=chromo, fontface="bold"),size=2.5, position = position_jitter(w = 0.3)) + guides(colour = "none")
    if(SSHFS) p_txt
    
    # outFile <- file.path(outFold, paste0(curr_var, "_matching_", tissue, "_consensus_boxplot_chromoDots.", plotType))
    # ggsave(plot=p_dot, file = outFile, width = widthBoxplot, height = heightBoxplot)
    # cat(paste0("... written: ", outFile, "\n"))
    # foo <- try(dev.off())
    
    outFile <- file.path(outFold, paste0(curr_var, "_matching_", tissue, "_consensus_boxplot_chromoLabs.", plotType))
    ggsave(plot=p_txt, file = outFile, width = widthBoxplot, height = heightBoxplot)
    cat(paste0("... written: ", outFile, "\n"))
    foo <- try(dev.off())
  }
  

} # end iterating var_to_plot

######################################################################################
cat(paste0("... written: ", logFile, "\n"))
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))






