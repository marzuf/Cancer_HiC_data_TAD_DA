startTime <- Sys.time()

SSHFS=F

# Rscript cmp_rebinning_custom_HiTC.R test_rao_matrix_10kb.txt test_rao_matrix_40kb.txt 10000 40000 chr1
# Rscript cmp_rebinning_custom_HiTC.R MCF-7/RAW_10kb/MC-7_chr16_RAW_10kb.hic.counts MCF-7/RAW_40kb/MC-7_chr16_RAW_40kb.hic.counts 10000 40000 chr16

# Rscript cmp_rebinning_custom_HiTC.R MCF-7/RAW_10kb/MC-7_chr16_RAW_10kb.hic.counts MCF-7/RAW_40kb/MC-7_chr16_RAW_40kb.hic.counts 10000 40000 chr16 > MCF-7/RAW_40kb/MC-7_chr16_RAW_40kb.hic.counts_checkRebinning.log

initBinFile= "test_rao_matrix_40kb.txt"
newBinFile="test_rao_matrix_40kb.txt"
oldBinSize=10000
newBinSize=40000
chromo="chr1"

args <- commandArgs(trailingOnly = TRUE)

stopifnot(length(args) == 5)

initBinFile <- args[1]
newBinFile <- args[2]

oldBinSize <- as.numeric(args[3])
newBinSize <- as.numeric(args[4])

chromo <- args[5]

stopifnot(file.exists(initBinFile))
stopifnot(file.exists(newBinFile))

logFile <- paste0(newBinFile, "_checkRebinning.log")
if(!SSHFS) file.remove(logFile)

stopifnot(!is.na(oldBinSize))
stopifnot(!is.na(newBinSize))

stopifnot(substr(chromo,1,3) == "chr")

suppressPackageStartupMessages(library(HiTC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(data.table, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(rtracklayer, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) # for preparing data for LGF
suppressPackageStartupMessages(library(Matrix, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) # for preparing data for LGF
# source("/mnt/etemp/marie/scripts/EZH2_old/EZH2_analysis/ezh2_utils_fct.R")
source("utils_fct.R")

##########################################################################################################################

# oldBinSize <- 10*10^3
# newBinSize <- 50*10^3

binInCoordinates <- TRUE
binfileHeader <- FALSE


######################## with HTC
htc_object <- createHTC(file=initBinFile, bin.size = oldBinSize, chr=chromo, dim = -1, reindex = binInCoordinates, header=binfileHeader)
hicMat_v0 <- as.matrix(intdata(htc_object))
dim(hicMat_v0)
htc_object <- binningC(htc_object, binsize=newBinSize, bin.adjust=FALSE, upa=TRUE, method="sum", optimize.by = "speed")
hicMat_v0b <- as.matrix(intdata(htc_object))
rownames(hicMat_v0b) <- colnames(hicMat_v0b) <- NULL

######################## with sparseMatrix
countDT <- fread(initBinFile, header=F, stringsAsFactors = F, col.names = c("binA", "binB", "count"))
countDT$idxA <- countDT$binA/oldBinSize + 1
countDT$idxB <- countDT$binB/oldBinSize + 1
hicMat_v1 <- sparseMatrix(i=countDT$idxA, j=countDT$idxB, x=countDT$count, symmetric=TRUE)
dim(hicMat_v1)
countDT$idxA <- countDT$idxB <- NULL

countDT$binA <- countDT$binA/oldBinSize
countDT$binB <- countDT$binB/oldBinSize

rebinDT <- rebin_sparseMatrix(sparseCountDT=countDT, initBinSize = oldBinSize, newBinSize=newBinSize, filled=FALSE)
stopifnot(all(colnames(rebinDT) == c("binA", "binB", "count")))

# rebinDT$idxA <- rebinDT$binA/newBinSize + 1
# rebinDT$idxB <- rebinDT$binB/newBinSize + 1
# head(rebinDT)
hicMat_v1b <- sparseMatrix(i=rebinDT$binA+1, j=rebinDT$binB+1, x=rebinDT$count, symmetric=TRUE)
dim(hicMat_v1b)


######################## COMPARISON HERE:
stopifnot(all.equal(as.matrix(hicMat_v0), as.matrix(hicMat_v1), check.attributes=F))
stopifnot(all.equal(as.matrix(hicMat_v0b), as.matrix(hicMat_v1b), check.attributes=F))

all.equal(as.matrix(hicMat_v0), as.matrix(hicMat_v1), check.attributes=F)
all.equal(as.matrix(hicMat_v0b), as.matrix(hicMat_v1b), check.attributes=F)

stopifnot(all.equal(as.matrix(hicMat_v0), as.matrix(hicMat_v1), check.attributes=F))
stopifnot(all.equal(as.matrix(hicMat_v0b), as.matrix(hicMat_v1b), check.attributes=F))


########################################################################################################### COMPARE WITH THE REBINNED DATA

countDT2 <- fread(newBinFile, header=F, stringsAsFactors = F, col.names = c("binA", "binB", "count"))

countDT2$idxA <- countDT2$binA/newBinSize + 1
countDT2$idxB <- countDT2$binB/newBinSize + 1
hicMat_v2 <- sparseMatrix(i=countDT2$idxA, j=countDT2$idxB, x=countDT2$count, symmetric=TRUE)

# hicMat_v2 <- symMatrix_from_sparseListDT(countDT2)
dim(hicMat_v2)

stopifnot(all.equal(as.matrix(hicMat_v0b), as.matrix(hicMat_v2), check.attributes=F))
stopifnot(all.equal(as.matrix(hicMat_v1b), as.matrix(hicMat_v2), check.attributes=F))

###
txt <- paste0("*** done ***\n")
printAndLog(txt, logFile)
txt <- paste0(startTime, "\n", Sys.time(), "\n")
printAndLog(txt, logFile)
cat(paste0("... written: ", logFile, "\n"))