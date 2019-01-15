library(data.table)
library(foreach)

options(scipen = 100)

startTime <- Sys.time()

printAndLog <- function(txt, logFile){
  cat(txt)
  cat(txt, file = logFile, append=T)
}

cat("> START check_matResol.R\n")
# Rscript check_matResol.R MCF-7

SSHFS <- FALSE

setDir <- ifelse(SSHFS, "/media/electron", "")
setDir <- ifelse(SSHFS, "~/media/electron", "")

args <- commandArgs(trailingOnly = TRUE)
if(length(args) == 1) {
  ds <- args[1]
}

binSizeKb <- 40
folderSuffix <- paste0("_", binSizeKb, "kb")
matrixSuffix <- "TopDom\\.matrix"

dsFold <- paste0(ds, folderSuffix)
stopifnot(dir.exists(dsFold))

outFold <- file.path("CHECK_MATRESOL", ds)
dir.create(outFold, recursive=TRUE)

logFile <- file.path(outFold, paste0(ds, "_check_matresol_logFile.txt"))
if(!SSHFS) file.remove(logFile)
if(SSHFS) logFile <- ""

txt <- paste0("!!! hard-coded binSizeKb = ", binSizeKb, "\n")
printAndLog(txt, logFile)

txt <- paste0("!!! hard-coded folderSuffix = ", folderSuffix, "\n")
printAndLog(txt, logFile)

txt <- paste0("!!! hard-coded matrixSuffix = ", matrixSuffix, "\n")
printAndLog(txt, logFile)

matFold <- file.path(dsFold, "TopDom_MAT")
stopifnot(dir.exists(matFold))

all_files <- list.files(matFold, full.names = TRUE, pattern = paste0(ds, "_chr.+", binSizeKb, "kb.+", matrixSuffix))
stopifnot(length(all_files) > 0)
stopifnot(length(all_files) <= 23)

chromoLevels <- paste0("chr", c(1:22, "X"))

# all_files=all_files[1:2]
check_resolDT <- foreach(inFile = all_files, .combine='rbind') %do% {
  
  txt <- paste0("> START file: ", basename(inFile), "\n")
  printAndLog(txt, logFile)
  
  inMat <- fread(inFile)
  
  coordDT <- inMat[,c(1:3)]
  
  stopifnot(length(unique(coordDT[,1])) == 1)
  stopifnot(length(unique(coordDT[,3]-coordDT[,2])) == 1)
  
  binSize <- unique(coordDT[,3]-coordDT[,2])
  stopifnot(binSize == (binSizeKb*1000))
  curr_chromo <- unique(coordDT[,1])
  
  txt <- paste0("... found chromo:\t", curr_chromo, "\n")
  printAndLog(txt, logFile)
  
  txt <- paste0("... found bin size:\t", binSize/1000, " kb", "\n")
  printAndLog(txt, logFile)
  
  hic_DT <- inMat[,-c(1:3)]
  stopifnot(nrow(hic_DT) == ncol(hic_DT))
  
  txt <- paste0("... ", curr_chromo, " - matrix dim.:\t", paste0(dim(hic_DT), collapse = " x "), "\n")
  printAndLog(txt, logFile)
  
  matrixRowSum <- rowSums(hic_DT, na.rm=T)
  
  outFile <- file.path(outFold, paste0(ds, "_", curr_chromo, "_matrixRowSum.Rdata"))
  save(matrixRowSum, file=outFile)
  cat(paste0("... written: ", outFile, "\n"))
  
  countSum <- sum(matrixRowSum)
  
  txt <- paste0("... matrix count sum:\t", round(countSum, 4), "\n")
  printAndLog(txt, logFile)
  
  txt <- "... summary matrix row sum:\n"
  printAndLog(txt, logFile)
  sink(logFile, append=T)
  print(summary(matrixRowSum))
  sink()
  
  rowAbove1000 <- sum(matrixRowSum >= 1000)/length(matrixRowSum)
  
  txt <- paste0("... # rows with >= 1000 counts:\t", sum(matrixRowSum >= 1000), "/", length(matrixRowSum), " (", round(rowAbove1000*100, 2),"%)\n")
  printAndLog(txt, logFile)
  txt <- "\n"
  printAndLog(txt, logFile)
  
  tmpDT <- data.frame(
    dataset = ds,
    chromo = curr_chromo,
    countSum = countSum,
    rowAbove1000 = rowAbove1000,
    stringsAsFactors = FALSE
  )
  colnames(tmpDT)[2] <- "chromo"
  
  tmpDT
}

check_resolDT$countSum <- round(check_resolDT$countSum, 4)
check_resolDT$rowAbove1000 <- round(check_resolDT$rowAbove1000, 4)

stopifnot(check_resolDT$chromo %in% chromoLevels)
check_resolDT$chromo <- factor(check_resolDT$chromo, levels = chromoLevels)
check_resolDT <- check_resolDT[order(as.numeric(check_resolDT$chromo), decreasing=F),]

outFile <- file.path(outFold, paste0(ds, "_check_resolDT.Rdata"))
save(check_resolDT, file=outFile)
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFold, paste0(ds, "_check_resolDT.txt"))
write.table(check_resolDT, col.names=T, row.names=F, sep="\t", quote=F, file = outFile)
cat(paste0("... written: ", outFile, "\n"))

######################################################################################
cat(paste0("... written: ", logFile, "\n"))
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))




