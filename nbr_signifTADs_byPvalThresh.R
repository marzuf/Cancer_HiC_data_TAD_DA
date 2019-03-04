
# Rscript nbr_signifTADs_byPvalThresh.R

startTime <- Sys.time()

cat("> START nbr_signifTADs_byPvalThresh.R \n")

require(foreach)
require(doMC)

SSHFS <- FALSE
setDir <- ifelse(SSHFS, "~/media/electron", "")

registerDoMC(ifelse(SSHFS, 2, 40))

outFold <- file.path("NBR_SIGNIFTADS_BYPVALTHRESH")
dir.create(outFold, recursive = TRUE)

plotCex <- 1.2
plotType <- "svg"
myHeight <- ifelse(plotType=="png", 400, 7)
myWidth <- ifelse(plotType=="png", 600, 10)

logFile=""

pval_step <- 0.01

pval_signif_thresh_seq <- seq(from=0, to=1, by=pval_step) 

#PIPELINE/OUTPUT_FOLDER/GSE105318_DLD1_40kb/TCGAcoad_msi_mss//emp_pval_combined.Rdata
script0_name <- "0_prepGeneData"
script11_name <- "11_runEmpPvalCombined"

pipoutFold <- file.path("PIPELINE", "OUTPUT_FOLDER")

all_hicexpr_ds <- unname(unlist(sapply(list.files(pipoutFold, full.names = TRUE), function(x) file.path(basename(x),  list.files(x)))))
stopifnot(dir.exists(file.path(pipoutFold, all_hicexpr_ds)))

cat("... start building signifTADs_allDS_data \n")
signifTADs_allDS_data <- foreach(ds = all_hicexpr_ds) %dopar% {
  cat("... start: ", ds, "\n")
  hicds <- dirname(ds)
  exprds <- basename(ds)
  stopifnot(dir.exists(hicds))
  dsPipOutDir <- file.path(pipoutFold, ds)
  stopifnot(dir.exists(dsPipOutDir))
    
  # RETRIEVE SIGNIF. TADs 
  stopifnot(dir.exists(file.path(dsPipOutDir, script11_name)))
  tad_pvalFile <- file.path(dsPipOutDir, script11_name, "emp_pval_combined.Rdata")
  stopifnot(file.exists(tad_pvalFile))
  tad_pvals <- eval(parse(text = load(tad_pvalFile)))
  adj_tad_pvals <- sort(p.adjust(tad_pvals, method="BH"))
  adj_tad_pvals
}  
  
names(signifTADs_allDS_data) <- all_hicexpr_ds

p_thresh <- 0.05
ratioSignif_by_thresh <- lapply(pval_signif_thresh_seq, function(p_thresh) {
  unlist(lapply(signifTADs_allDS_data, function(x) mean(x <= p_thresh) ))
})
names(ratioSignif_by_thresh) <- pval_signif_thresh_seq
  
nSignif_by_thresh <- lapply(pval_signif_thresh_seq, function(p_thresh) {
  unlist(lapply(signifTADs_allDS_data, function(x) sum(x <= p_thresh) ))
})
names(nSignif_by_thresh) <- pval_signif_thresh_seq

nSignif_by_thresh_DT <- as.data.frame(nSignif_by_thresh)
ratioSignif_by_thresh_DT <- as.data.frame(ratioSignif_by_thresh)

myTit <- paste0("# signif. TADs vs. p-val. thresh.")
mySub <- paste0("(pval_step = ", pval_step, ")")

### plot 1) # signif TADs
myxlab <- paste0("nbr signif. TADs")

outFile <- file.path(outFold, paste0("nbr_signif_by_pval.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
boxplot(nSignif_by_thresh_DT, las=2,
        cex.lab = plotCex, cex.axis = plotCex,
        xlab = myxlab,
        main = myTit)
mtext(text = mySub, side = 3)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

### plot 2) # signif TADs (log10)
myxlab <- paste0("nbr signif. TADs (log10)")

outFile <- file.path(outFold, paste0("nbr_signif_by_pval_log10.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
boxplot(log10(nSignif_by_thresh_DT), las=2,
        cex.lab = plotCex, cex.axis = plotCex,
        xlab = myxlab,
        main = myTit)
mtext(text = mySub, side = 3)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

### plot 3) ratio signif TADs
myxlab <- paste0("ratio signif. TADs")
outFile <- file.path(outFold, paste0("ratio_signif_by_pval.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
boxplot(ratioSignif_by_thresh_DT, las=2,
        xlab = myxlab,
        cex.lab = plotCex, cex.axis = plotCex,
        main = myTit)
mtext(text = mySub, side = 3)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

