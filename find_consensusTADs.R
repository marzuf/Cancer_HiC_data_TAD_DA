options(scipen = 100)

library(doMC)
library(foreach)

startTime <- Sys.time()

cat("> START find_consensusTADs.R\n")
# Rscript find_consensusTADs.R  MCF-7_40kb ENCSR549MGQ_T47D_40kb

SSHFS <- FALSE

registerDoMC(ifelse(SSHFS, 2, 40))

setDir <- ifelse(SSHFS, "/media/electron", "")

if(SSHFS) setwd("/media/electron/mnt/etemp/marie/Dixon2018_integrative_data")

source("utils_fct.R")

source("../Dixon2018_integrative_data/wsm_TADconsensus_withCoverage_version2_fct_withAllTxt.R")


#*****************************

# as input: should be able to retrieve _final_domains.txt
# automatically try to retrieve as much chromos as possible
# MCF-7_40kb/FINAL_DOMAINS/MCF-7_chr2_KR_40kb_final_domains.txt 
# 
# otuput should look like
# <cl1cl2>_40kb/FINAL_DOMAINS/<cl1cl2>_chr15_KR_40kb_final_domains.txt

# command line: all main input folders
# Rscript find_consensusTADs.R cl1_40kb cl2_40kb cl3_40kb ...

all_inFolders <- c("MCF-7_40kb", "Caki2_40kb")
all_inFolders <- c("MCF-7_40kb", "ENCSR549MGQ_T47D_40kb")


all_inFolders <- commandArgs(trailingOnly = TRUE)
stopifnot(file.exists(all_inFolders))


folderSuffix <- unique(gsub(".+(_.+?kb)$", "\\1", all_inFolders))
stopifnot(length(folderSuffix) == 1)
all_cl <- gsub(folderSuffix, "", all_inFolders)
consensusClName <- paste0(all_cl, collapse="")
outFolder <- file.path(paste0(consensusClName, folderSuffix), "FINAL_DOMAINS")
stopifnot(length(outFolder) == 1)

dir.create(outFolder, recursive=TRUE)


#*****************************
logFile <- file.path(outFolder, paste0(consensusClName, "_consensusTADs.logFile.txt"))
if(!SSHFS) system(paste0("rm -f ", logFile))
if(SSHFS) logFile <- ""

# hard-coded settings for the consensus
txt  <- paste0("*** HARD-CODED PIPELINE SETTINGS:\n")
printAndLog(txt, logFile)
domainFolder <- "FINAL_DOMAINS"
txt  <- paste0("... domainFolder\t=\t", domainFolder, "\n")
tadfile_suffixPattern <- "_final_domains.txt"
txt  <- paste0("... tadfile_suffixPattern\t=\t", tadfile_suffixPattern, "\n")

# hard-coded settings for the consensus
txt  <- paste0("*** HARD-CODED SETTINGS FOR CONSENSUS TAD CALLING ALGORITHM:\n")
printAndLog(txt, logFile)

set_tolRad <- 80000
set_conservThresh <- ifelse(length(all_inFolders) == 3, 2/3, 
                            ifelse(length(all_inFolders) == 2, 1, NA))
if(is.na(set_conservThresh)) stop("error: cannot infer set_conservThresh\n")
set_weightValue <- NULL
set_coverageThresh <-  1
set_tadfileHeader <-  FALSE
set_chrSize <- NULL
set_ncpu <- ifelse(SSHFS, 2, 40)
set_fileAsInput <- TRUE
set_logFile <- NULL

txt  <- paste0("... set_tolRad\t=\t", set_tolRad, "\n")
printAndLog(txt, logFile)
txt  <- paste0("... set_conservThresh\t=\t", set_conservThresh, "\n")
printAndLog(txt, logFile)
txt  <- paste0("... set_weightValue\t=\t", set_weightValue, "\n")
printAndLog(txt, logFile)
txt  <- paste0("... set_coverageThresh\t=\t", set_coverageThresh, "\n")
printAndLog(txt, logFile)
txt  <- paste0("... set_tadfileHeader\t=\t", as.character(set_tadfileHeader), "\n")
printAndLog(txt, logFile)
txt  <- paste0("... set_chrSize\t=\t", as.character(set_chrSize), "\n")
printAndLog(txt, logFile)
txt  <- paste0("... set_ncpu\t=\t", set_ncpu, "\n")
printAndLog(txt, logFile)
txt  <- paste0("... set_fileAsInput\t=\t", set_fileAsInput, "\n")
printAndLog(txt, logFile)
txt  <- paste0("... set_logFile\t=\t", as.character(set_logFile), "\n")
printAndLog(txt, logFile)

txt  <- paste0("*** SETTINGS RETRIEVED FROM COMMAND LINE:\n")
printAndLog(txt, logFile)

txt  <- paste0("... all_inFolders\t=\t", paste0(all_inFolders, collapse=", "), "\n")
printAndLog(txt, logFile)

###############################################################################################################
###############################################################################################################
###############################################################################################################

domain_inFolders <- file.path(all_inFolders, domainFolder)
stopifnot(dir.exists(domain_inFolders))

# retrieve for each CL the available domain files (hold in list to retrieve common chromo)
domain_inFiles <- lapply(domain_inFolders, function(inFolder) {
  inFiles <- list.files(inFolder, full.names = TRUE, recursive=FALSE, pattern=paste0(tadfile_suffixPattern, "$"))
  stopifnot(length(inFiles) > 0)
  inFiles
})
stopifnot(length(domain_inFiles) == length(all_cl))

# retrieve normalization (for name of output file)
normMeth <- unique(sapply( unlist(domain_inFiles), function(dfile) gsub(paste0(".+_chr.+?_(.+?)", folderSuffix, tadfile_suffixPattern), "\\1", basename(dfile))))
stopifnot(length(normMeth) == 1)

# retrieve for each CL the available chromos
all_chromos <- lapply(seq_along(domain_inFiles), function(i) {
  cl_chromos <- unique(gsub(paste0(".+_(chr.+?)_.+", folderSuffix, tadfile_suffixPattern), "\\1", basename(domain_inFiles[[i]])))
  stopifnot(cl_chromos %in% c(paste0("chr", 1:23), "chrX"))
  cl_chromos
})

# retrieve common chromo
intersectChromo <- Reduce(intersect, all_chromos)
stopifnot(length(intersectChromo) > 0)
stopifnot(intersectChromo %in% c(paste0("chr", 1:23), "chrX"))

txt <- paste0("... found common chromo(s):\t", paste0(intersectChromo, collapse=","), "\n")
printAndLog(txt, logFile)

# intersectChromo <- intersectChromo[1]
# intersectChromo <- "chrX"

for(chromo in intersectChromo) {
  
  cat(paste0("> START chromo = ", chromo, "\n"))
  
  chromo_files <- lapply(seq_along(domain_inFiles), function(i) {
    curr_files <- domain_inFiles[[i]]
    curr_chr_file <- curr_files[grep(paste0(chromo, "_"), basename(curr_files))]
    stopifnot(file.exists(curr_chr_file))
    stopifnot(length(curr_chr_file) == 1)
    curr_chr_file
  })
  
    domainsDT <- get_consensus_TADs(chromo_files,
                                    tolRad=set_tolRad,
                                    conservThresh=set_conservThresh,
                                    weightValue = set_weightValue,
                                    coverageThresh = set_coverageThresh,
                                  tadfileHeader = set_tadfileHeader,
                                  chrSize=set_chrSize,
                                  ncpu=set_ncpu,
                                  fileAsInput=set_fileAsInput,
                                  logFile=set_logFile)
  
    
  cat("nrow consensus:", nrow(domainsDT), "\n")
  
  if(nrow(domainsDT) > 1 ){
    for(i in 2:nrow(domainsDT)) {
      # cat(paste0(i, "\t"))
      stopifnot( domainsDT[i, 2] > domainsDT[i-1, 3]  )
    }
  }
  # <cl1cl2>_40kb/FINAL_DOMAINS/<cl1cl2>_chr15_KR_40kb_final_domains.txt
  # outFile <- file.path(outFolder, paste0(consensusClName, "_", chromo,"_", normMeth, "_", folderSuffix, "_", tadfile_suffixPattern))
  outFile <- file.path(outFolder, paste0(consensusClName, "_", chromo,"_", normMeth,  folderSuffix, tadfile_suffixPattern))
  write.table(domainsDT, file = outFile, sep="\t", quote=F, col.names=F, row.names=F)
  cat(paste0("... written:\t", outFile, "\n"))  
}
 

######################################################################################
cat(paste0("... written: ", logFile, "\n"))
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))




