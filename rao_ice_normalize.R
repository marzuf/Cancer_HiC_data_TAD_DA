startTime <- Sys.time()

SSHFS=F

# Rscript cmp_rebinning_custom_HiTC.R \
#GSE73782_PC3/RAW_40kb/GSE73782_PC3_chr21_RAW_10kb.hic.counts \
#GSE73782_PC3_40kb/NORM_MAT_ICE/GSE73782_PC3_chr21_ICE_40kb.hic.matrix chr20 40000

# Rscript rao_ice_normalize.R GSE73782_PC3/RAW_40kb/GSE73782_PC3_chr21_RAW_40kb.hic.counts GSE73782_PC3_40kb/NORM_MAT_ICE/GSE73782_PC3_chr21_ICE_40kb.hic.matrix chr21 40000

args <- commandArgs(trailingOnly = TRUE)

stopifnot(length(args) == 4)

inFile <- args[1]
outFile <- args[2]
chromo <- args[3]
binSize <- as.numeric(args[4])

stopifnot(file.exists(inFile))

dir.create(dirname(outFile), recursive = TRUE)

logFile <- file.path(dirname(outFile), paste0(basename(outFile), "_ICE_normalize.log"))
if(!SSHFS) file.remove(logFile)

stopifnot(!is.na(binSize))

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

ICE_maxiter <- 1000

######################## load map with HTC
htc_object <- createHTC(file=inFile, bin.size = binSize, 
                        chr=chromo, dim = -1, reindex = binInCoordinates, 
                        header=binfileHeader)

htc_object_ICE <- normICE(htc_object, max_iter=ICE_maxiter)

htc_ice_matrix <- as.data.frame(as.matrix(intdata(htc_object_ICE)))

colsTopDom <- data.frame(chromo=chromo, 
                         binA=seq(0,length.out = nrow(htc_ice_matrix), by=binSize), 
                         binB=seq(binSize, length.out = nrow(htc_ice_matrix), by = binSize),
                         stringsAsFactors = FALSE)

outMat <- cbind(colsTopDom, as.data.frame(as.matrix(htc_ice_matrix)))

cat("... write outFile\n")
write.table(outMat, file = outFile, col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
cat("... written: ",  outFile, "\n")


###
txt <- paste0("*** done ***\n")
printAndLog(txt, logFile)
txt <- paste0(startTime, "\n", Sys.time(), "\n")
printAndLog(txt, logFile)
cat(paste0("... written: ", logFile, "\n"))