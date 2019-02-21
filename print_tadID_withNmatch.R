# Rscript intersect_topTADs_acrossDS.R 3

startTime <- Sys.time()

cat("> START print_tadID_withNmatch.R \n")

# Rscript print_tadID_withNmatch.R 21
# Rscript print_tadID_withNmatch.R 21 13 INTERSECT_topTADs_ACROSSDS/top3/all_matchDT.Rdata

SSHFS <- FALSE
setDir <- ifelse(SSHFS, "~/media/electron", "")

plotType <- "svg"

require(ggplot2)
source( file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg"), "set_dataset_colors.R"))
head(score_DT)
source( file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_11_18"), "analysis_utils.R"))
source( file.path("colors_utils.R"))
dataset_proc_colors <- setNames(score_DT$proc_col, score_DT$dataset)
length(dataset_proc_colors)

source("utils_fct.R")

outFolder <- file.path("PRINT_TADID_WITHNMATCH")
dir.create(outFolder)

nMatch=21
all_nMatch <- c(21, 13)

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) > 0)
all_nMatch <- as.numeric(args)
if(is.na(args[length(args)])) {
  matchDT_file <- all_nMatch[length(all_nMatch)]
}else{
  matchDT_file <- file.path("INTERSECT_topTADs_ACROSSDS", "top3", "all_matchDT.Rdata")
}

if(SSHFS) {
  checkFile <- ""
  
} else {
  checkFile <- file.path(outFolder, paste0("matching_with_nMatch_", paste(args, collapse="_"), ".txt"))
  file.remove(checkFile)
  
}

stopifnot(file.exists(matchDT_file))
all_matchDT <- eval(parse(text = load(matchDT_file)))

all_bestMatchDT <- do.call(rbind,
                           lapply(split(all_matchDT,list(all_matchDT$query_id,all_matchDT$matching_hicds,all_matchDT$matching_exprds ),drop=T), 
                                  function(subDT) subDT[which.max(subDT$matchingRatio),]))
rownames(all_bestMatchDT) <- NULL

all_nMatchDT <- aggregate(matching_id ~ query_id, FUN=length, data = all_bestMatchDT)
head(all_nMatchDT)
colnames(all_nMatchDT) <- c("query_id", "all_nMatch")

# queryID_matchDT <- merge(all_nMatchDT, diffExprds_nMatchDT, by="query_id", all.x = TRUE, all.y=TRUE)
queryID_matchDT <- all_nMatchDT
stopifnot(!is.na(queryID_matchDT$all_nMatch))
# sum(is.na(queryID_matchDT$diffExprds_nMatch))
stopifnot(!duplicated(queryID_matchDT$query_id))

srtd_DT <- queryID_matchDT[order(queryID_matchDT$all_nMatch, decreasing = T),]

nMatch = 21
for(nMatch in all_nMatch){
  mytxt <- paste0("#*********************************************** For nMatch = ", nMatch, "\n")
  stopifnot(!is.na(nMatch))
  printAndLog(txt = mytxt, logFile = checkFile)
  
  mytxt <- paste0("\n#***> TAD IDs:", "\n")
  printAndLog(txt = mytxt, logFile = checkFile)
  nMatch_DT <- srtd_DT[srtd_DT$all_nMatch == nMatch,]
  write.table(nMatch_DT[, "query_id", drop=F], file="", sep="\t", quote=F, col.names=F, row.names = FALSE)
  write.table(nMatch_DT[, "query_id", drop=F], file=checkFile, append=TRUE, sep="\t", quote=F, col.names=F, row.names = FALSE)
  
  queryID_nmatch <- nMatch_DT$query_id
  queryID_matching <- sapply(queryID_nmatch, function(x) {
    paste0(sort(c(x, all_bestMatchDT$matching_id[all_bestMatchDT$query_id==x])), collapse=",\n")
  } )
  
  mytxt <- paste0("\n#***> unique set(s):", "\n")
  printAndLog(txt = mytxt, logFile = checkFile)
  unique_sets <- unique(queryID_matching)
  mytxt <- paste0("> ", unique_sets, sep="\n")
  printAndLog(txt = mytxt, logFile = checkFile)
  
  mytxt <- paste0("\n\n#***> TAD IDs with matching TADs:")#, "\n")
  printAndLog(txt = mytxt, logFile = checkFile)
  
  mytxt <- paste(paste0("\n> ", names(queryID_matching), "\t:"), queryID_matching, sep="\n")
  printAndLog(txt = mytxt, logFile = checkFile)
  
  txt <- paste0("\n#", rep("#", 50), "\n")
  printAndLog(txt = mytxt, logFile = checkFile)

}
cat(paste0("... written: ", checkFile, "\n"))


# ######################################################################################
# ######################################################################################
# ######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))
