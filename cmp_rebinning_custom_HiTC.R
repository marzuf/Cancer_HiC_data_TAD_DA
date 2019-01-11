startTime <- Sys.time()

### CHECKED:
# Rscript check_rebinning.R ALL_29_11_10kb/KARPAS/DMSO/COUNTS/KARPAS_DMSO_chr1_10000_aggregCounts.txt ALL_29_11_10kb/KARPAS/DMSO/COUNTS/KARPAS_DMSO_chr1_10000_aggregCounts.txt
# (see list of checked datasets in check_rebinning.txt !)
args <- commandArgs(trailingOnly = TRUE)

stopifnot(length(args) == 2)


# initBinFile <- "../ALL_29_11_10kb/KARPAS/DMSO/COUNTS/KARPAS_DMSO_chr1_10000_aggregCounts.txt"
# newBinFile <- "../ALL_29_11_10kb/KARPAS/DMSO/COUNTS/KARPAS_DMSO_chr1_50000_aggregCounts.txt"
initBinFile <- args[1]
newBinFile <- args[2]

chromo <- gsub(".+_.+_(.+?)_.+", "\\1", basename(newBinFile))

stopifnot(substr(chromo,1,3) == "chr")

suppressPackageStartupMessages(library(HiTC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(data.table, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(rtracklayer, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) # for preparing data for LGF
suppressPackageStartupMessages(library(Matrix, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) # for preparing data for LGF
source("../EZH2_analysis/ezh2_utils_fct.R")

##########################################################################################################################

# binSize <- 10*10^3
# newBinSize <- 50*10^3

binSize <- as.numeric(gsub(".+_(.+?)_aggregCounts.txt", "\\1", basename(initBinFile)))
stopifnot(!is.na(binSize))
newBinSize <- as.numeric(gsub(".+_(.+?)_aggregCounts.txt", "\\1", basename(newBinFile)))
stopifnot(!is.na(newBinSize))

binInCoordinates <- FALSE
binfileHeader <- FALSE

resolThresh_percent <- 0.8
resolThresh_count <- 1000

######################## with HTC
htc_object <- createHTC(file=initBinFile, bin.size = binSize, chr=chromo, dim = -1, reindex = binInCoordinates, header=binfileHeader)
hicMat_v0 <- as.matrix(intdata(htc_object))
dim(hicMat_v0)
htc_object <- binningC(htc_object, binsize=newBinSize, bin.adjust=FALSE, upa=TRUE, method="sum", optimize.by = "speed")
hicMat_v0b <- as.matrix(intdata(htc_object))
rownames(hicMat_v0b) <- colnames(hicMat_v0b) <- NULL

######################## with sparseMatrix
countDT <- fread(initBinFile, header=F, stringsAsFactors = F, col.names = c("bin1", "bin2", "count"))
hicMat_v1 <- symMatrix_from_sparseListDT(countDT)
dim(hicMat_v1)
rebinDT <- rebin_sparseMatrix(sparseCountDT=countDT, initBinSize = binSize, newBinSize=newBinSize, filled=FALSE)
stopifnot(all(colnames(rebinDT) == c("bin1", "bin2", "count")))
hicMat_v1b <- symMatrix_from_sparseListDT(rebinDT) 

######################## COMPARISON HERE:
stopifnot(all.equal(as.matrix(hicMat_v0), as.matrix(hicMat_v1), check.attributes=F))
stopifnot(all.equal(as.matrix(hicMat_v0b), as.matrix(hicMat_v1b), check.attributes=F))

all.equal(as.matrix(hicMat_v0), as.matrix(hicMat_v1), check.attributes=F)
all.equal(as.matrix(hicMat_v0b), as.matrix(hicMat_v1b), check.attributes=F)

stopifnot(all.equal(as.matrix(hicMat_v0), as.matrix(hicMat_v1), check.attributes=F))
stopifnot(all.equal(as.matrix(hicMat_v0b), as.matrix(hicMat_v1b), check.attributes=F))

checkCriterion(hicMat = hicMat_v0b, percentThresh = resolThresh_percent, countThresh = resolThresh_count, logF=NULL)
checkCriterion(hicMat = hicMat_v1b, percentThresh = resolThresh_percent, countThresh = resolThresh_count, logF=NULL)


########################################################################################################### COMPARE WITH THE REBINNED DATA

countDT2 <- fread(newBinFile, header=F, stringsAsFactors = F, col.names = c("bin1", "bin2", "count"))
hicMat_v2 <- symMatrix_from_sparseListDT(countDT2)
dim(hicMat_v2)

stopifnot(all.equal(as.matrix(hicMat_v0b), as.matrix(hicMat_v2), check.attributes=F))
stopifnot(all.equal(as.matrix(hicMat_v1b), as.matrix(hicMat_v2), check.attributes=F))

###
cat("*** done ***\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))
