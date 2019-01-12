#!/usr/bin/env Rscript

startTime <- Sys.time()
cat(paste0("*** START prepFromHiC_prepFile.R - ", startTime, "\n"))

suppressPackageStartupMessages(library(optparse, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(data.table, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

options(scipen=100)


SSHFS  <- F
setDir <- ifelse(SSHFS, "~/media/electron", "")

source("utils_fct.R")

# Rscript prepFromHiC_prepFile.R -f countFile -F preFile -c chromo -b binSize
# Rscript prepFromHiC_prepFile.R -f MCF-7/RAW_40kb/MC-7_chr21_RAW_40kb.hic.counts -F MCF-7/PRE_40kb/MC-7_chr21_RAW_40kb.pre -c chr21 -b binSize

# countFile="MCF-7/RAW_10kb/MC-7_chr16_RAW_10kb.hic.counts"
# chromo = 16
# binSize = 40000

option_list = list(
  
  make_option(c("-c", "--chromo"), type="character", default=NULL,
              help="chromosome", metavar="character"),
  
  make_option(c("-b", "--binSize"), type="integer", default=NULL,
              help="bin size", metavar="integer"),
  
  make_option(c("-f", "--countFile"), type="character", default=NULL,
              help="path to input file", metavar="character"),
  
  make_option(c("-F", "--preFile"), type="character", default=NULL,
              help="path to output file", metavar="character")

)

opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);
if(
  is.null(opt$countFile) | is.null(opt$preFile)|
  is.null(opt$chromo) | is.null(opt$binSize)){
  stop("Error - missing input argument !\n")
}

chromo <- opt$chromo
countFile <- opt$countFile
preFile <- opt$preFile
binSize <- opt$binSize

stopifnot(file.exists(countFile))

logFile <- paste0(preFile, "_preFile.log")
if(!SSHFS) file.remove(logFile)

dir.create(dirname(preFile), recursive = TRUE)

cat("countFile = ", countFile, "\n")
cat("preFile = ", preFile, "\n")
cat("chromo = ", chromo, "\n")
cat("binSize = ", binSize, "\n")


cat("... read countFile\n")
curr_dt <- fread(countFile, header=FALSE, col.names=c("binA", "binB", "count")) 
stopifnot(ncol(curr_dt) == 3)
head(curr_dt)
stopifnot(is.numeric(unlist(curr_dt$count)))
stopifnot(is.numeric(unlist(curr_dt$binA)))
stopifnot(is.numeric(unlist(curr_dt$binB)))

### prepare the datable for the pre file
fragVec1 <- 0
fragVec2 <- 1
posVec1 <- curr_dt$binA
posVec2 <- curr_dt$binB
countVec <- curr_dt$count
chrVec <- chromo
preDT <- data.frame(str1=fragVec1, chr1=chrVec, pos1=posVec1, frag1=fragVec1, 
                    str1=fragVec1, chr1=chrVec, pos1=posVec2, frag1=fragVec2, counts=countVec, stringsAsFactors = FALSE)

cat(paste0("... write \"pre\" to file \n"))
write.table(preDT, file=preFile, row.names = F, col.names=F, sep=" ", quote=F)
cat(paste0("... written: ", preFile, "\n"))

chrSize <- max(c(curr_dt$binA, curr_dt$binB), na.rm=T)
# is 0-based
chrSize <- chrSize+binSize

# need to write chromosome size file

chrSizeDT <- data.frame(
  chromo = chromo,
  chrSize=chrSize,
  stringsAsFactors = FALSE
)
outFile <- file.path(dirname(preFile), paste0(chromo, ".size"))
write.table(chrSizeDT, file=outFile, row.names = F, col.names=F, sep="\t", quote=F)
cat(paste0("... written: ", outFile, "\n"))

txt <- paste0("*** done ***\n")
printAndLog(txt, logFile)
txt <- paste0(startTime, "\n", Sys.time(), "\n")
printAndLog(txt, logFile)
cat(paste0("... written: ", logFile, "\n"))


    
    
    # cat(paste0("......... run juicer tools pre \n"))
    # # add -n to not normalize
    # # -r binSize -> to do only for a single resolution
    # # -c chromo -> only for the chromo
    # # -d -> only intrachromo
    # command <- paste("java -Xmx2g -jar", path_to_juicer_tools, "pre -n -d -r", binSize, "-c", chromo, countFile_pre, preFile_hic, countFile_size)
    # cat(paste0(command, "\n"))
    # system(command)
    # 

