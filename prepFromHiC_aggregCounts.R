#!/usr/bin/env Rscript

startTime <- Sys.time()
cat(paste0("*** START prepFromHiC_aggregCounts.R - ", startTime, "\n"))

suppressPackageStartupMessages(library(optparse, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(data.table, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

options(scipen=100)

############################################################################# 
############################################################################# HARD-CODED SETTINGS
############################################################################# 
SSHFS  <- F
setDir <- ifelse(SSHFS, "~/media/electron", "")
source(paste0("utils_fct.R"))
########################################################################################################################################################## 

# Rscript prepFromHiC_aggregCounts.R -f inFile -F outFile -b 10000 -B 40000




option_list = list(
  
   
  ### !!!!! INTERCHROMO !!!!!!!!!!!! or not ???
  ### !!!!! CHECK THAT IT SHOULD BE ZERO-BASED
  
  
  make_option(c("-c", "--chromo"), type="character", default=NULL,
              help="chromosome", metavar="character"),
  
  
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
if( is.null(opt$chromo) | 
   is.null(opt$inFile) | is.null(opt$outFile)|
   is.null(opt$inBinSize) |   is.null(opt$outBinSize) ){
  stop("Error - missing input argument !\n")
}

oldBinSize <- opt$inBinSize
newBinSize <- opt$outBinSize
inFile <- opt$inFile
outFile <- opt$outFile
chromo <- opt$chromo

stopifnot(file.exists(inFile))


stopifnot(newBinSize %% 1 == 0)
stopifnot(binSize %% 1 == 0)
stopifnot(newBinSize >= binSize)
stopifnot((newBinSize/binSize) %% 1 == 0)


dir.create(dirname(outFile), recursive = TRUE)


txt <- paste0("*** START ", chromo, "\n")
cat(txt)

cat("inFile = ", mainDir, "\n")
cat("outFile = ", outFile, "\n")
cat("oldBinSize = ", oldBinSize, "\n")
cat("newBinSize = ", newBinSize, "\n")


cat("... read inFile\n")
curr_dt <- fread(inFile, header=FALSE, col.names=c("binA", "binB", "counts")) 
stopifnot(ncol(curr_dt) == 3)
head(curr_dt)
stopifnot(is.numeric(unlist(curr_dt$counts)))
stopifnot(is.numeric(unlist(curr_dt$binA)))
stopifnot(is.numeric(unlist(curr_dt$binB)))
# ensure all binsB > binsA !
stopifnot(all(curr_dt$binB >= curr_dt$binA))


# convert to the same format as for the bin files from Stephanie
# 0 -> 0, 50000 -> 1
# hard-coded check: data from GSE63525 are 50 kb!

curr_dt$binA <- curr_dt$binA/binSize
curr_dt$binB <- curr_dt$binB/binSize
stopifnot(all(curr_dt$binA %%1 == 0))
stopifnot(all(curr_dt$binB %%1 == 0))

if(oldBinSize == newBinSize) {
  cat("... no rebinning needed\n") 
  rebinDT <- curr_dt
} else{
  cat(paste0("... start rebinning from ", oldBinSize, " to ", newBinSize, "\n"))
  colnames(curr_dt) <- c("bin1", "bin2", "count")
  rebinDT <- rebin_sparseMatrix(sparseCountDT=curr_dt, 
                                initBinSize=oldBinSize, 
                                newBinSize=newBinSize, 
                                filled=FALSE) 
}

# GSE63525: to get the same format for all GSE63525 data -> multiply by the binSize
cat("... convert back coordinates and write new matrix to file\n")
# write with new coordinates
rebinDT$bin1 <- rebinDT$bin1 * newBinSize
rebinDT$bin2 <- rebinDT$bin2 * newBinSize

write.table(rebinDT, file = outFile, quote=F, sep="\t", col.names=F, row.names=F)

cat(paste0("> FINISHED for ", chromo, " \n... written: ", outFile, "\n"))
    

cat("*** DONE \n")
cat(paste0("... check logFile: ", logFile, "\n"))
cat(paste0(startTime, "\n", Sys.time(), "\n"))


