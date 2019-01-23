
# Rscript print_topTADs.R <hicds> <exprds> <nTop>

# Rscript print_topTADs.R GSE105318_DLD1_40kb TCGAcoad_msi_mss 20
# Rscript print_topTADs.R GSE105318_DLD1_40kb TCGAcoad_msi_mss 20

# Rscript print_topTADs.R ENCSR444WCZ_A549_40kb TCGAluad_mutKRAS_mutEGFR 20
# Rscript print_topTADs.R NCI-H460_40kb TCGAluad_mutKRAS_mutEGFR 20
# Rscript print_topTADs.R ENCSR444WCZ_A549NCI-H460_40kb TCGAluad_mutKRAS_mutEGFR 20
# Rscript print_topTADs.R TCGAluad_mutKRAS_mutEGFR 20

# Rscript print_topTADs.R ENCSR079VIJ_G401_40kb TCGAkich_norm_kich 20
# Rscript print_topTADs.R ENCSR401TBQ_Caki2_40kb TCGAkich_norm_kich 20
# Rscript print_topTADs.R ENCSR079VIJ_G401ENCSR401TBQ_Caki2_40kb TCGAkich_norm_kich 20
# Rscript print_topTADs.R TCGAkich_norm_kich 20

# Rscript print_topTADs.R ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_wt_mutCTNNB1 20
# Rscript print_topTADs.R ENCSR862OGI_RPMI-7951_40kb TCGAskcm_wt_mutCTNNB1 20
# Rscript print_topTADs.R ENCSR312KHQ_SK-MEL-5ENCSR862OGI_RPMI-7951_40kb TCGAskcm_wt_mutCTNNB1 20
# Rscript print_topTADs.R TCGAskcm_wt_mutCTNNB1 20

args <- commandArgs(trailingOnly = TRUE)

hicds <- "GSE105318_DLD1_40kb"
exprds <- "TCGAcoad_msi_mss"
topThresh <- 5
if(length(args) == 3){
  hicds <- args[1]
  exprds <- args[2]
  nTop <- as.numeric(args[3])
  
} else if(length(args) == 2){
  hicds <- ""
  exprds <- args[1]
  nTop <- as.numeric(args[2])
  
} else {
  stop("... invalid argument #\n")
}

stopifnot(!is.na(nTop))
stopifnot(nTop > 0)

#PIPELINE/OUTPUT_FOLDER/GSE105318_DLD1_40kb/TCGAcoad_msi_mss//emp_pval_combined.Rdata
script11_name <- "11_runEmpPvalCombined"

if(nchar(hicds) > 0) {
  stopifnot(dir.exists(hicds))
  pipOutDir <- file.path("PIPELINE", "OUTPUT_FOLDER", hicds, exprds)
} else {
  pipOutDir <- file.path("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom", "OUTPUT_FOLDER", hicds, exprds)
}

stopifnot(dir.exists(pipOutDir))
stopifnot(dir.exists(file.path(pipOutDir, script11_name)))

tad_pvalFile <- file.path(pipOutDir, script11_name, "emp_pval_combined.Rdata")
stopifnot(file.exists(tad_pvalFile))
tad_pvals <- eval(parse(text = load(tad_pvalFile)))

adj_tad_pvals <- sort(p.adjust(tad_pvals, method="BH"))


if(topThresh > 1) {
  topThresh <- min(c(topThresh, length(adj_tad_pvals)))
  pvalThresh <- as.numeric(adj_tad_pvals[topThresh])
  stopifnot(!is.na(pvalThresh))
  stopifnot(pvalThresh <= 1)
} else {
    pvalThresh <- topThresh
    topThresh <- min(which(adj_tad_pvals <= pvalThresh))
}

top_pvals <- adj_tad_pvals[adj_tad_pvals <= topThresh]
stopifnot(!is.na(top_pvals))
top_tads <- names(top_pvals)
top_pvals <- as.numeric(top_pvals)

resultDT <- data.frame(
  # top_pvals = round(top_pvals, 4),
  # top_pvals = top_pvals,
  top_pvals = sprintf("%.2e",top_pvals),
  top_tads = top_tads,
  stringsAsFactors = FALSE
)

cat(paste0("********** ", hicds, " - ", exprds, "\n\n"))
cat(paste0("> topThresh\t=\t", topThresh, "\n"))
cat(paste0("> pvalThresh\t=\t", sprintf("%.2e",pvalThresh), "\n\n"))
cat(paste0("*** Top ranking TADs (adj. emp. p-val. combined): \n"))

#write.table(head(resultDT), col.names=TRUE, row.names=FALSE, sep="\t", quote=F, file = "")

write.table(resultDT, col.names=TRUE, row.names=FALSE, sep="\t", quote=F, file = "")









