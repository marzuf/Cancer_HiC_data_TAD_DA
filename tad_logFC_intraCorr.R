
# Rscript tad_logFC_intraCorr.R # default 0.05
# Rscript tad_logFC_intraCorr.R 0.05


startTime <- Sys.time()

cat("> START tad_logFC_intraCorr.R \n")

SSHFS <- FALSE

require(foreach)
require(doMC)

source("utils_fct.R")

plotType <- "svg"
myHeight <- ifelse(plotType=="png", 500, 7)
myWidth <- myHeight

printVar <- function(x){
  cat(paste0(x, " = ", eval(parse(text=x)), "\n"))
}
writeSepToFile <- function(mysymbol, myfile, ...){
  printAndLog(txt = paste0("\n", paste0(rep(mysymbol, 50), collapse=""), "\n"), logFile = checkFile)
}

writeValueToFile <- function(myvarname, myfile, ...) {
  txt <- paste0(myvarname, " =\t", eval(parse(text=myvarname)), "\n")
  cat(txt, file="")
  cat(txt, file = myfile, append = TRUE, ...)
}

writeTableToFile <- function(mytable, myfile) {
  write.table(mytable, col.names=TRUE, row.names=FALSE, append=T, sep="\t", quote=F, file=myfile)
  write.table(mytable, col.names=TRUE, row.names=FALSE, append=T, sep="\t", quote=F, file="")
}

registerDoMC(ifelse(SSHFS, 2, 40))

build_signifTADs_allDS_data <- TRUE


setDir <- ifelse(SSHFS, "~/media/electron", "")

signifThresh <- 0.05


args <- commandArgs(trailingOnly = TRUE)
if(length(args) > 0) {
signifThresh <- as.numeric(args[1])
stopifnot(!is.na(signifThresh))
}

stopifnot(signifThresh >= 0)
stopifnot(signifThresh <= 1)

outFolder <- file.path("TAD_LOGFC_INTRACORR", paste0("signif", signifThresh))
dir.create(outFolder, recursive = TRUE)

# checkFile <- file.path(outFolder, paste0("check_file_traceback_", toptad_id, "_", bottad_id, ".txt"))
# file.remove(checkFile)
checkFile <- file.path(outFolder, paste0("check_file_traceback", ".txt"))
file.remove(checkFile)


logFile=""

stopifnot(!is.na(signifThresh))
stopifnot(signifThresh > 0)

#PIPELINE/OUTPUT_FOLDER/GSE105318_DLD1_40kb/TCGAcoad_msi_mss//emp_pval_combined.Rdata
script0_name <- "0_prepGeneData"
script9_name <- "9_runEmpPvalMeanTADLogFC"
script10_name <- "10_runEmpPvalMeanTADCorr"
script11_name <- "11_runEmpPvalCombined"

tad_match_folder <- file.path("INTERSECT_topTADs_ACROSSDS", "top3")
stopifnot(dir.exists(tad_match_folder))

all_bestMatchDT_file <- file.path(tad_match_folder, "all_bestMatchDT.Rdata")
stopifnot(file.exists(all_bestMatchDT_file))

signifTADs_allDS_data_file <- file.path(tad_match_folder, "signifTADs_allDS_data.Rdata")
stopifnot(file.exists(signifTADs_allDS_data_file))

all_matchDT_file <- file.path(tad_match_folder, "all_matchDT.Rdata")
stopifnot(file.exists(all_matchDT_file))

ratio_matchingSignifTAD_DT_file <- file.path(tad_match_folder, "ratio_matchingSignifTAD_DT.Rdata")
stopifnot(file.exists(ratio_matchingSignifTAD_DT_file))

hicds_exprds_asMatch_DT_file <- file.path(tad_match_folder, "hicds_exprds_asMatch_DT.Rdata")
stopifnot(file.exists(hicds_exprds_asMatch_DT_file))


pipOutFolder <- file.path("PIPELINE", "OUTPUT_FOLDER")

all_hicexpr_ds <- unname(unlist(sapply(list.files(pipOutFolder, full.names = TRUE), function(x) file.path(basename(x),  list.files(x)))))
stopifnot(dir.exists(file.path(pipOutFolder, all_hicexpr_ds)))

ds=all_hicexpr_ds[1]

if(build_signifTADs_allDS_data){
  
cat("... start building allPvals_allDS_DT data \n")
  
allPvals_allDS_DT <- foreach(ds = all_hicexpr_ds, .combine='rbind') %dopar% {
  
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
    tad=all_tads,
    tad_id=paste(hicds, exprds, all_tads, sep="_"),
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
# load("TAD_LOGFC_INTRACORR/signif0.05/allPvals_allDS_DT.Rdata")
# hicds             exprds adj_pvalFC adj_pvalCorr adj_pvalComb
# 1 ENCSR079VIJ_G401_40kb TCGAkich_norm_kich 0.02913459  0.002099790 2.915007e-05
# 2 ENCSR079VIJ_G401_40kb TCGAkich_norm_kich 0.02913459  0.002950338 3.456350e-05

xvar <- "adj_pvalFC"
yvar <- "adj_pvalCorr"
myTit <- "all TCGA datasets - all TADs"
mySub <- paste0("(# TADs = ", length(allPvals_allDS_DT$tad_id) , ")")

plotCex <- 1.2

outFile <- file.path(outFolder, paste0("pvalCorr_vs_pvalFC_allTADs", ".", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
myx <- allPvals_allDS_DT[,xvar]
myy <- allPvals_allDS_DT[,yvar]
densplot(x = myx,
         xlab=paste0(xvar),
         xlim=range(myx,na.rm=T),
         ylim=range(myy,na.rm=T),
         y = myy,
         ylab=paste0(yvar),
         main = myTit,
         cex=0.7,
         cex.axis=plotCex,
         cex.lab=plotCex,
         pch=16
)
mtext(side=3, text=mySub)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFolder, paste0("pvalCorr_vs_pvalFC_allTADs_log10", ".", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
myx <- allPvals_allDS_DT[,xvar]
myy <- allPvals_allDS_DT[,yvar]
densplot(x = myx,
         xlab=paste0(xvar),
         xlim=range(myx,na.rm=T),
         ylim=range(myy,na.rm=T),
         y = myy,
         ylab=paste0(yvar),
         main = myTit,
         cex=0.7,
         cex.axis=plotCex,
         cex.lab=plotCex,
         pch=16
)
mtext(side=3, text=mySub)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


signifPvals_allDS_DT <- allPvals_allDS_DT[allPvals_allDS_DT$adj_pvalFC <= signifThresh & 
                                            allPvals_allDS_DT$adj_pvalCorr <= signifThresh,]

stopifnot(!duplicated(signifPvals_allDS_DT$tad_id))

myTit <- paste0("all TCGA datasets - signif. TADs (pval <= ", signifThresh, ")")
mySub <- paste0("(# TADs = ", length(signifPvals_allDS_DT$tad_id) , ")")

myx <- signifPvals_allDS_DT[,xvar]
myy <- signifPvals_allDS_DT[,yvar]

outFile <- file.path(outFolder, paste0("pvalCorr_vs_pvalFC_signifTADs", ".", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
densplot(x = myx,
         xlab=paste0(xvar),
         xlim=c(min(myx,na.rm=T),signifThresh),
         ylim=c(min(myy,na.rm=T),signifThresh),
         y = myy,
         ylab=paste0(yvar),
         main = myTit,
         cex=0.7,
         cex.axis=plotCex,
         cex.lab=plotCex,
         pch=16
)
mtext(side=3, text=mySub)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

myx <- -log10(signifPvals_allDS_DT[,xvar])
myy <-  -log10(signifPvals_allDS_DT[,yvar])

outFile <- file.path(outFolder, paste0("pvalCorr_vs_pvalFC_signifTADs_log10", ".", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
densplot(x = myx,
         xlab=paste0(xvar),
         xlim=c(min(myx,na.rm=T),signifThresh),
         ylim=c(min(myy,na.rm=T),signifThresh),
         y = myy,
         ylab=paste0(yvar),
         main = myTit,
         cex=0.7,
         cex.axis=plotCex,
         cex.lab=plotCex,
         pch=16
)
mtext(side=3, text=mySub)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))




# ###################################################################################### ADD the main TAD
load("TAD_LOGFC_INTRACORR/signif0.05/allPvals_allDS_DT.Rdata")
head(allPvals_allDS_DT)

all_bestMatchDT <- eval(parse(text = load(all_bestMatchDT_file)))
head(all_bestMatchDT)

all_nMatchDT <- aggregate(matching_id ~ query_id, FUN=length, data = all_bestMatchDT)
head(all_nMatchDT)
colnames(all_nMatchDT) <- c("query_id", "all_nMatch")

srtd_DT <- all_nMatchDT[order(all_nMatchDT$all_nMatch, decreasing = TRUE),]

top_to_plot <- 5

xvar_all <- "adj_pvalFC"
yvar_all <- "adj_pvalCorr"

xvar_top <- "adj_pvalFC"
yvar_top <- "adj_pvalCorr"


# datasetwise
# outFile <- file.path(outFolder, paste0("pvalCorr_vs_pvalFC_signifTADs", ".", plotType))
# do.call(plotType, list(outFile, height=myHeight, width=myWidth))
# plot(x = allPvals_allDS_DT[,xvar],
#          xlab=paste0(xvar),
#          xlim=c(0,1),
#          ylim=c(0,1),
#          y = allPvals_allDS_DT[,yvar],
#          ylab=paste0(yvar),
#          main = myTit,
#          cex=plotCex/2,
#          cex.axis=plotCex,
#          cex.lab=plotCex,
#      col="grey",
#          pch=16
# )

i=1
all_x_to_plot <- c()
all_y_to_plot <- c()
all_labels_to_plot <- c()
all_cols_to_plot <- c()
all_nmatch_legend <- c()

plotted_sets <- list()

nPlotted <- i<-  0
while(nPlotted < top_to_plot) {
  
  i <- i+1
  
  cat("... nPlotted = ", nPlotted, "\n")
  
  plot_list <- list()
  
  tad_id <- srtd_DT$query_id[i]
  nmatch_legend <- srtd_DT$all_nMatch[i]
  
  matching_ids <- all_bestMatchDT$matching_id[all_bestMatchDT$query_id == tad_id]
  stopifnot(!duplicated(matching_ids))
  
  to_plot_set <- c(tad_id, matching_ids)
  stopifnot(!duplicated(to_plot_set))
  
  if(any(unlist(lapply(plotted_sets, function(x) all(to_plot_set %in% x))))) {
    next
  } else {
    plotted_sets <- c(plotted_sets, list(to_plot_set))
    stopifnot(list(to_plot_set) %in% plotted_sets)
    nPlotted <- nPlotted+1
  }
  
  stopifnot(any(unlist(lapply(plotted_sets, function(x) length(setdiff(to_plot_set, x)) == 0))))
  
  tmpDT <- allPvals_allDS_DT
  stopifnot(!duplicated(tmpDT$tad_id))
  rownames(tmpDT) <- tmpDT$tad_id
  stopifnot(to_plot_set %in% rownames(tmpDT) )
  
  to_plot_x <- tmpDT[to_plot_set, xvar]
  to_plot_y <- tmpDT[to_plot_set, yvar]
  to_plot_labels <- paste(tmpDT[to_plot_set, "hicds"],tmpDT[to_plot_set, "exprds"], tmpDT[to_plot_set, "tad"], sep="\n")
  to_plot_col <- length(plotted_sets)
  all_x_to_plot <- c(all_x_to_plot, to_plot_x)
  all_y_to_plot <- c(all_y_to_plot, to_plot_y)
  all_labels_to_plot <- c(all_labels_to_plot, to_plot_labels)
  all_cols_to_plot <- c(all_cols_to_plot, to_plot_col)
  all_nmatch_legend <- c(all_nmatch_legend, nmatch_legend)
  
  
  #datasetwise
  # points(x=to_plot_x, y = to_plot_y, col =to_plot_col , pch=4)
  # text(x=to_plot_x,y=to_plot_y,col=to_plot_col, labels = to_plot_labels, cex = 0.7)
  
}




outFile <- file.path(outFolder, paste0("pvalCorr_vs_pvalFC_matchTADs_", "nTop",top_to_plot , ".", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))

plot(x = allPvals_allDS_DT[,xvar],
     xlab=paste0(xvar),
     xlim=range(all_x_to_plot),
     ylim=range(all_y_to_plot),
     y = allPvals_allDS_DT[,yvar],
     ylab=paste0(yvar),
     main = myTit,
     cex=plotCex/2,
     cex.axis=plotCex,
     cex.lab=plotCex,
     col="grey",
     pch=16
)
points(x=all_x_to_plot, y = all_y_to_plot, col =all_cols_to_plot , pch=4)
text(x=all_x_to_plot,y=all_y_to_plot,col=all_cols_to_plot, labels = all_labels_to_plot, cex = 0.7)
legend("topleft", legend = all_nmatch_legend, col = all_cols_to_plot, pch=4, bty="n")
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


outFile <- file.path(outFolder, paste0("pvalCorr_vs_pvalFC_matchTADs_", "nTop",top_to_plot , "_log10.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))

all_x_to_plot_log10 <- -log10(all_x_to_plot)
all_y_to_plot_log10 <- -log10(all_y_to_plot)

plot(x =  -log10(allPvals_allDS_DT[,xvar]),
     xlab=paste0(xvar, "[-log10]"),
     xlim=range(all_x_to_plot_log10),
     ylim=range(all_y_to_plot_log10),
     y =  -log10(allPvals_allDS_DT[,yvar]),
     ylab=paste0(yvar, "[-log10]"),
     main = myTit,
     cex=plotCex/2,
     cex.axis=plotCex,
     cex.lab=plotCex,
     col="grey",
     pch=16
)
points(x=all_x_to_plot_log10, y = all_y_to_plot_log10, col =all_cols_to_plot , pch=4)
text(x=all_x_to_plot_log10,y=all_y_to_plot_log10,col=all_cols_to_plot, labels = all_labels_to_plot, cex = 0.7)
legend("topleft", legend = all_nmatch_legend, col = all_cols_to_plot, pch=4, bty="n")
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

# ###################################################################################### check the algo
# load("TAD_LOGFC_INTRACORR/signif0.05/allPvals_allDS_DT.Rdata")
# head(allPvals_allDS_DT)

# saved from tad_logFC_intraCorr.R 
# save(signifTADs_allDS_data, file = outFile)
# save(all_matchDT, file = outFile)
# save(all_bestMatchDT, file = outFile)
# save(ratio_matchingSignifTAD_DT, file = outFile)
# save(hicds_exprds_asMatch_DT, file = outFile)

# 1) check table1
signifTADs_allDS_data <- eval(parse(text = load(signifTADs_allDS_data_file)))
head(signifTADs_allDS_data)

# 2) check table2
all_matchDT <- eval(parse(text = load(all_matchDT_file)))
head(all_matchDT)

# 3) check table3
all_bestMatchDT <- eval(parse(text = load(all_bestMatchDT_file)))
head(all_bestMatchDT)

# 4) check table4
ratio_matchingSignifTAD_DT <- eval(parse(text = load(ratio_matchingSignifTAD_DT_file)))
head(ratio_matchingSignifTAD_DT)

# 5) check table5
hicds_exprds_asMatch_DT <- eval(parse(text = load(hicds_exprds_asMatch_DT_file)))
head(hicds_exprds_asMatch_DT)

# traceback 2 TADs (from top and bottom to check the script !)
toptad_exprds <- "TCGAskcm_lowInf_highInf"
toptad_hicds <- "ENCSR312KHQ_SK-MEL-5_40kb"
toptad_tad <-"chr6_TAD71"
toptad_id <- paste(toptad_hicds, toptad_exprds, toptad_tad, sep="_")


bottad_exprds <- "TCGAlusc_norm_lusc"
bottad_hicds <- "NCI-H460_40kb-MEL-5_40kb"
bottad_tad <-"chr1_TAD232"
bottad_id <- paste(bottad_hicds, bottad_exprds, bottad_tad, sep="_")



#### START CHECK HERE
curr_exprds=toptad_exprds
curr_hicds=toptad_hicds
curr_tad=toptad_tad
curr_id=toptad_id

all_exprds=c(toptad_exprds, bottad_exprds)
all_hicds=c(toptad_hicds, bottad_hicds)
all_tad=c(toptad_tad, bottad_tad)
all_id=c(toptad_id, bottad_id)


for(i in 1:2) {
  signifTADs_allDS_data <- eval(parse(text = load(signifTADs_allDS_data_file)))
  
  
  
  curr_exprds=all_exprds[i]
  curr_hicds=all_hicds[i]
  curr_tad=all_tad[i]
  curr_id=all_id[i]
  # Error: curr_id %in% curr_signifTADs_allDS_data[["ids"]] is not TRUE
  
  writeSepToFile(mysymbol="*", myfile = checkFile)
  
  writeValueToFile(myvarname="curr_exprds", myfile=checkFile)
  writeValueToFile(myvarname="curr_hicds", myfile=checkFile)
  writeValueToFile(myvarname="curr_tad", myfile=checkFile)
  writeValueToFile(myvarname="curr_id", myfile=checkFile)
  
  cat(file.path(curr_hicds, curr_exprds))
  cat("\n")
  cat(names(signifTADs_allDS_data))
  
  x1=file.path(curr_hicds, curr_exprds) 
  save(x1, file="x1.Rdata")
  x2=names(signifTADs_allDS_data)
save(x2, file="x2.Rdata")
load("x1.Rdata")
load("x2.Rdata")

  stopifnot(file.path(curr_hicds, curr_exprds) %in% names(signifTADs_allDS_data))

  curr_signifTADs_allDS_data <- signifTADs_allDS_data[[file.path(curr_hicds, curr_exprds)]]
  
  cat(length(curr_signifTADs_allDS_data))
  cat("\n")
  
  # "ids"     cat("\n")
   # "idDT"     "geneDT"   "signifDT" "posDT"    "matchDT" 
  cat(head(curr_signifTADs_allDS_data[["ids"]]))
  cat("\n")
  
  cat(curr_id)
  cat("\n")
  head(curr_signifTADs_allDS_data[["ids"]])
  head(names(curr_signifTADs_allDS_data))
  
  stopifnot(curr_id %in% curr_signifTADs_allDS_data[["ids"]])
  # curr_ids <- curr_signifTADs_allDS_data[["ids"]][curr_signifTADs_allDS_data[["ids"]] == curr_id]
  # stopifnot(length(curr_ids) > 0)
  
  head(curr_signifTADs_allDS_data[["idDT"]])
  stopifnot(curr_id %in% curr_signifTADs_allDS_data[["idDT"]]$ID)
  curr_idDT <- curr_signifTADs_allDS_data[["idDT"]][curr_signifTADs_allDS_data[["idDT"]]$ID == curr_id,]
  stopifnot(nrow(curr_idDT) > 0)
  
  head(curr_signifTADs_allDS_data[["geneDT"]])
  stopifnot(curr_id %in% curr_signifTADs_allDS_data[["geneDT"]]$ID)
  curr_geneDT <- curr_signifTADs_allDS_data[["geneDT"]][curr_signifTADs_allDS_data[["geneDT"]]$ID == curr_id,]
  stopifnot(nrow(curr_geneDT) > 0)
  
  head(curr_signifTADs_allDS_data[["signifDT"]])
  stopifnot(curr_id %in% curr_signifTADs_allDS_data[["signifDT"]]$ID)
  curr_signifDT <- curr_signifTADs_allDS_data[["signifDT"]][curr_signifTADs_allDS_data[["signifDT"]]$ID == curr_id,]
  stopifnot(nrow(curr_signifDT) > 0)
  
  head(curr_signifTADs_allDS_data[["posDT"]])
  stopifnot(curr_id %in% curr_signifTADs_allDS_data[["posDT"]]$ID)
  curr_posDT <- curr_signifTADs_allDS_data[["posDT"]][curr_signifTADs_allDS_data[["posDT"]]$ID == curr_id,]
  stopifnot(nrow(curr_posDT) > 0)
  
  head(curr_signifTADs_allDS_data[["matchDT"]])
  stopifnot(curr_id %in% curr_signifTADs_allDS_data[["matchDT"]]$ID)
  curr_matchDT <- curr_signifTADs_allDS_data[["matchDT"]][curr_signifTADs_allDS_data[["matchDT"]]$ID == curr_id,]
  stopifnot(nrow(curr_matchDT) > 0)
  
  
  writeSepToFile(mysymbol="*", myfile = checkFile)
  
  writeTableToFile(mytable=curr_idDT, myfile=checkFile)
  writeSepToFile(mysymbol="*", myfile = checkFile)
  
  writeTableToFile(mytable=curr_geneDT, myfile=checkFile)
  writeSepToFile(mysymbol="*", myfile = checkFile)
  
  writeTableToFile(mytable=curr_signifDT, myfile=checkFile)
  writeSepToFile(mysymbol="*", myfile = checkFile)
  
  writeTableToFile(mytable=curr_posDT, myfile=checkFile)
  writeSepToFile(mysymbol="*", myfile = checkFile)
  
  writeTableToFile(mytable=curr_matchDT, myfile=checkFile)
}
cat(paste0("... written: ", checkFile, "\n"))

# ######################################################################################
# ######################################################################################
# ######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))
