
# Rscript tad_logFC_intraCorr.R 0.05

startTime <- Sys.time()

cat("> START tad_logFC_intraCorr.R \n")

SSHFS <- FALSE

require(foreach)
require(doMC)

source("utils_fct.R")

printVar <- function(x){
  cat(paste0(x, " = ", eval(parse(text=x)), "\n"))
}

registerDoMC(ifelse(SSHFS, 2, 40))

build_signifTADs_allDS_data <- TRUE


setDir <- ifelse(SSHFS, "~/media/electron", "")

signifThresh <- 0.05


args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) > 0)
signifThresh <- as.numeric(args[1])

if(signifThresh == 1) {
  warning("signifThresh == 1 is ambiguous; will be considered as a nTop not pval thresh !\n")
}

outFolder <- file.path("TAD_LOGFC_INTRACORR", paste0("signif", signifThresh))
dir.create(outFolder, recursive = TRUE)

logFile=""

stopifnot(!is.na(signifThresh))
stopifnot(signifThresh > 0)

#PIPELINE/OUTPUT_FOLDER/GSE105318_DLD1_40kb/TCGAcoad_msi_mss//emp_pval_combined.Rdata
script0_name <- "0_prepGeneData"
script9_name <- "9_runEmpPvalMeanTADLogFC"
script10_name <- "10_runEmpPvalMeanTADCorr"
script11_name <- "11_runEmpPvalCombined"

pipOutFolder <- file.path("PIPELINE", "OUTPUT_FOLDER")



all_hicexpr_ds <- unname(unlist(sapply(list.files(pipOutFolder, full.names = TRUE), function(x) file.path(basename(x),  list.files(x)))))
stopifnot(dir.exists(file.path(pipOutFolder, all_hicexpr_ds)))



ds=all_hicexpr_ds[1]

if(build_signifTADs_allDS_data){
  
cat("... start building allPvals_allDS_DT data \n")
  
allPvals_allDS_DT <- foreach(ds = all_hicexpr_ds) %dopar% {
  
  cat("... start: ", ds, "\n")
  
  hicds <- dirname(ds)
  exprds <- basename(ds)
  stopifnot(dir.exists(hicds))
  dsPipOutDir <- file.path(pipOutFolder, ds)
  stopifnot(dir.exists(dsPipOutDir))
  
  # RETRIEVE logFC pval
  stopifnot(dir.exists(file.path(dsPipOutDir, script9_name)))
  tad_pvalFCFile <- file.path(dsPipOutDir, script9_name, "emp_pval_meanLogFC.Rdata")
  stopifnot(file.exists(tad_pvalFCFile))
  tad_pvalFC <- eval(parse(text = load(tad_pvalFCFile))) # not adjusted
  tad_pvalFC <- sort(p.adjust(tad_pvalFC, method="BH"))
  
  # tad_pvalFC <- tad_pvalFC[tad_pvalFC <= signifThresh]
  
  # RETRIEVE TADcorr pval
  stopifnot(dir.exists(file.path(dsPipOutDir, script10_name)))
  tad_pvalCorrFile <- file.path(dsPipOutDir, script10_name, "emp_pval_meanCorr.Rdata")
  stopifnot(file.exists(tad_pvalCorrFile))
  tad_pvalCorr <- eval(parse(text = load(tad_pvalCorrFile)))
  tad_pvalCorr <- sort(p.adjust(tad_pvalCorr, method="BH"))
  
  # tad_pvalCorr <- tad_pvalCorr[tad_pvalCorr <= signifThresh]
  
  # RETRIEVE SIGNIF. TADs 
  stopifnot(dir.exists(file.path(dsPipOutDir, script11_name)))
  
  tad_pvalFile <- file.path(dsPipOutDir, script11_name, "emp_pval_combined.Rdata")
  stopifnot(file.exists(tad_pvalFile))
  tad_pvals <- eval(parse(text = load(tad_pvalFile)))
  tad_pvalComb <- sort(p.adjust(tad_pvals, method="BH"))
  
  # tad_pvalComb <- tad_pvalComb[tad_pvalComb <= signifThresh]
  
  all_tads <- Reduce(intersect, list(names(tad_pvalFC), names(tad_pvalCorr), names(tad_pvalComb)))
  stopifnot(length(all_tads) == length(tad_pvalFC))
  stopifnot(length(all_tads) == length(tad_pvalCorr))
  stopifnot(length(all_tads) == length(tad_pvalComb))
  
  outDT <- data.frame(
    hicds=hicds,
    exprds=exprds,
    adj_pvalFC = as.numeric(tad_pvalFC[all_tads]),
    adj_pvalCorr = as.numeric(tad_pvalCorr[all_tads]),
    adj_pvalComb = as.numeric(tad_pvalComb[all_tads]),
    stringsAsFactors = FALSE
  )
  stopifnot(!is.na(outDT))
  outDT
}
outFile <- file.path(outFolder, "allPvals_allDS_DT.Rdata")
save(allPvals_allDS_DT, file = outFile)
cat(paste0("... written: ", outFile, "\n"))

} else { # if not to build
  outFile <- file.path(outFolder, "allPvals_allDS_DT.Rdata")
  allPvals_allDS_DT <- eval(parse(text = load(outFile)))
}
  
head(allPvals_allDS_DT)
