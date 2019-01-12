#!/usr/bin/env Rscript

startTime <- Sys.time()
cat(paste0("*** START prepFromHiC_aggregCounts.R - ", startTime, "\n"))

suppressPackageStartupMessages(library(optparse, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(data.table, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(Matrix, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))


options(scipen=100)

############################################################################# 
############################################################################# HARD-CODED SETTINGS
############################################################################# 
SSHFS  <- F
setDir <- ifelse(SSHFS, "~/media/electron", "")
source(paste0("utils_fct.R"))
########################################################################################################################################################## 

# Rscript prepFromHiC_aggregCounts.R -f inFile -F outFile -b 10000 -B 40000

#  Rscript prepFromHiC_aggregCounts.R -f MCF-7/RAW_10kb/MC-7_chr16_RAW_10kb.hic.counts -F MCF-7/RAW_40kb/MC-7_chr16_RAW_40kb.hic.counts -b 10000 -B 40000
#  Rscript prepFromHiC_aggregCounts.R -f MCF-7/RAW_10kb/MC-7_chr16_RAW_10kb.hic.counts -F MCF-7/RAW_40kb/MC-7_chr16_RAW_40kb.hic.counts -b 10000 -B 40000 > MCF-7/RAW_40kb/MC-7_chr16_RAW_40kb.hic.counts_aggregCounts.log

#  Rscript prepFromHiC_aggregCounts.R -f test_rao_matrix_10kb.txt -F test_rao_matrix_40kb.txt -b 10000 -B 40000


# inFile="MCF-7/RAW_10kb/MC-7_chr16_RAW_10kb.hic.counts"
# outFile="MCF-7/RAW_40kb/MC-7_chr16_RAW_40kb.hic.counts"

inFile="test_rao_matrix_10kb.txt"
newBinSize=40000
oldBinSize=10000

option_list = list(
  
   
  ### !!!!! INTERCHROMO !!!!!!!!!!!! or not ???
  ### !!!!! CHECK THAT IT SHOULD BE ZERO-BASED
  
  
  
  make_option(c("-f", "--inFile"), type="character", default=NULL,
              help="path to input file", metavar="character"),
  
  make_option(c("-F", "--outFile"), type="character", default=NULL,
              help="path to output file", metavar="character"),
  
  
  make_option(c("-b", "--inBinSize"), type="integer", default=NULL,
              help="init size of the bins", metavar="integer"),
  
  make_option(c("-B", "--outBinSize"), type="integer", default=NULL,
              help="new size of the bins", metavar="integer")

  
)

opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);
if(
   is.null(opt$inFile) | is.null(opt$outFile)|
   is.null(opt$inBinSize) |   is.null(opt$outBinSize) ){
  stop("Error - missing input argument !\n")
}

oldBinSize <- opt$inBinSize
newBinSize <- opt$outBinSize
inFile <- opt$inFile
outFile <- opt$outFile

stopifnot(file.exists(inFile))


logFile <- paste0(inFile, "_aggregCounts.log")
if(!SSHFS) file.remove(logFile)


stopifnot(newBinSize %% 1 == 0)
stopifnot(oldBinSize %% 1 == 0)
stopifnot(newBinSize >= oldBinSize)
stopifnot((newBinSize/oldBinSize) %% 1 == 0)


dir.create(dirname(outFile), recursive = TRUE)

cat("inFile = ", inFile, "\n")
cat("outFile = ", outFile, "\n")
cat("oldBinSize = ", oldBinSize, "\n")
cat("newBinSize = ", newBinSize, "\n")


cat("... read inFile\n")
curr_dt <- fread(inFile, header=FALSE, col.names=c("binA", "binB", "count")) 
stopifnot(ncol(curr_dt) == 3)
head(curr_dt)
stopifnot(is.numeric(unlist(curr_dt$count)))
stopifnot(is.numeric(unlist(curr_dt$binA)))
stopifnot(is.numeric(unlist(curr_dt$binB)))
# ensure all binsB > binsA !
stopifnot(all(curr_dt$binB >= curr_dt$binA))

# convert to the same format as for the bin files from Stephanie
# 0 -> 0, 50000 -> 1
# hard-coded check: data from GSE63525 are 50 kb!

curr_dt$binA <- curr_dt$binA/oldBinSize
curr_dt$binB <- curr_dt$binB/oldBinSize
stopifnot(all(curr_dt$binA %%1 == 0))
stopifnot(all(curr_dt$binB %%1 == 0))

if(oldBinSize == newBinSize) {
  cat("... no rebinning needed\n") 
  rebinDT <- curr_dt
} else{
  cat(paste0("... start rebinning from ", oldBinSize, " to ", newBinSize, "\n"))
  rebinDT <- rebin_sparseMatrix(sparseCountDT=curr_dt, 
                                initBinSize=oldBinSize, 
                                newBinSize=newBinSize, 
                                filled=FALSE) 
}


inDT <- as.data.frame(curr_dt)
inMat <- as.matrix(sparseMatrix(i=(inDT$binA+1), j=(inDT$binB+1), x=inDT$count, symmetric=TRUE))
# write.table(inMat, file = "inMat.txt", quote=F, sep="\t", col.names=F, row.names=F)

outDT <- as.data.frame(rebinDT)
outMat <- as.matrix(sparseMatrix(i=(outDT$binA+1), j=(outDT$binB+1), x=outDT$count, symmetric=TRUE))
# write.table(outMat, file = "outMat.txt", quote=F, sep="\t", col.names=F, row.names=F)

stopifnot(sum(inMat) == sum(outMat))


# GSE63525: to get the same format for all GSE63525 data -> multiply by the binSize
cat("... convert back coordinates and write new matrix to file\n")
# write with new coordinates
rebinDT$binA <- rebinDT$binA * newBinSize
rebinDT$binB <- rebinDT$binB * newBinSize

write.table(rebinDT, file = outFile, quote=F, sep="\t", col.names=F, row.names=F)

cat(paste0("... written: ", outFile, "\n"))
    
txt <- paste0("*** done ***\n")
printAndLog(txt, logFile)
txt <- paste0(startTime, "\n", Sys.time(), "\n")
printAndLog(txt, logFile)
cat(paste0("... written: ", logFile, "\n"))

