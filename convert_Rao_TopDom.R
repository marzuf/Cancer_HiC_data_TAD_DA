startTime <- Sys.time()

script_name <- "convert_Rao_TopDom"

cat(paste0("> START ", script_name, ".R\n"))

suppressPackageStartupMessages(library(data.table, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(Matrix, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

options(scipen=100)

SSHFS <- F
setDir <- ifelse(SSHFS, "~/media/electron", "")

if(SSHFS) setwd(file.path(setDir, "/mnt/etemp/marie/Cancer_HiC_data_TAD_DA"))

source("utils_fct.R")

# Rscript convert_Rao_TopDom.R inFile outFile chromo binSizeBp
# inFile format:
  # 1st col=bin 1 in kb (0-based)
  # 2st col=bin 2 in kb (0-based)
  # 3d col=normalized count
    # 100000  100000  1487.574
    # 750000  750000  1050.8768
    # 800000  800000  799.052

# inFile="PC3_rep_12_KR_50kb.hic.matrix"
# binSizeBp=50000
# chromo="chr1"
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 4)
inFile <- args[1]
outFile <- args[2]
chromo <- args[3]
binSizeBp <- as.numeric(args[4])
stopifnot(file.exists(inFile))

dir.create(dirname(outFile), recursive = TRUE)

logFile <- file.path(dirname(outFile), paste0(script_name, "_logFile.txt"))
if(!SSHFS & file.exists(logFile)) file.remove(logFile)

txt <- paste0("*** RETRIEVED FROM COMMAND LINE:\n")
printAndLog(txt, logFile)
txt <- paste0("...> inFile = ", inFile, "\n")
printAndLog(txt, logFile)
txt <- paste0("...> outFile = ", outFile, "\n")
printAndLog(txt, logFile)
txt <- paste0("...> binSizeBp = ", binSizeBp, "\n")
printAndLog(txt, logFile)

cat("... read inFile\n")
inDT <- fread(inFile, col.names=c("binA", "binB", "count"))
stopifnot(is.numeric(inDT$binA))
stopifnot(is.numeric(inDT$binB))
stopifnot(is.numeric(inDT$count))
stopifnot(inDT$binA%%binSizeBp == 0)
stopifnot(inDT$binB%%binSizeBp == 0)

cat("... remove NA count (would make TopDom crash) ... \n")
nbrNA <- sum(is.na(inDT$count))
initNrow <- nrow(inDT)
inDT <- na.omit(inDT)
newNrow <- nrow(inDT)
# txt <- paste0("...> initNrow = ", initNrow, "\n")
# printAndLog(txt, logFile)
# txt <- paste0("...> newNrow = ", newNrow, "\n")
# printAndLog(txt, logFile)
# txt <- paste0("...> nbrNA = ", nbrNA, "\n")
# printAndLog(txt, logFile)
stopifnot((initNrow-newNrow) == nbrNA)
txt <- paste0("... # NA removed: ", nbrNA, "\n")
printAndLog(txt, logFile)

inDT$idxA <- inDT$binA/binSizeBp + 1
inDT$idxB <- inDT$binB/binSizeBp + 1

# check because will be done symmetric

# integer vectors of the same length specifying the locations (row and column indices) of the non-zero (or non-TRUE) entries of the matrix
# inMat <- sparseMatrix(i=c(1,1,1,2,2,2,3,3,3), j=c(1,2,3,1,2,3,1,2,3), x = c(1:9)) # 1-based !
# inMat <- sparseMatrix(i=c(1,1,1,2,2,2,3,3,3), j=c(1,2,3,1,2,3,1,2,3), x = c(1:9)) # 1-based !
# inMat <- sparseMatrix(i=c(1,2,2,3,3,3), j=c(1,1,2,1,2,3), x = c(1:6), symmetric = TRUE) # 1-based !
# inMat

cat("... build symmetric matrix\n")
inMat <- sparseMatrix(i=inDT$idxA, j=inDT$idxB, x=inDT$count, symmetric=TRUE)

stopifnot(isSymmetric(inMat))

txt <- paste0("... found matrix dim: ", paste0(dim(inMat), collapse="x"), "\n")
printAndLog(txt, logFile)

colsTopDom <- data.frame(chromo=chromo, 
                         binA=seq(0,length.out = nrow(inMat), by=binSizeBp), 
                         binB=seq(binSizeBp, length.out = nrow(inMat), by = binSizeBp),
                         stringsAsFactors = FALSE)

outMat <- cbind(colsTopDom, as.data.frame(as.matrix(inMat)))

cat("... write outFile\n")
write.table(outMat, file = outFile, col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
cat("... written: ",  outFile, "\n")


######################################################################################
######################################################################################
######################################################################################
cat("*** DONE\n")
cat("... written: ",  logFile, "\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))


