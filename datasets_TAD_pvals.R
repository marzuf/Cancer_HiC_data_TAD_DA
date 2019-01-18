startTime <- Sys.time()
cat(paste0("> START datasets_TAD_pvals.R\n"))

# Rscript datasets_TAD_pvals.R

source("utils_fct.R")
options(scipen=100)

suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
suppressPackageStartupMessages(library(flux, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 


SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

computeAUC <- TRUE
buildTable <- TRUE

plotType <- "svg"
myHeightGG <- 7
myWidthGG <- 10
myHeight <- ifelse(plotType == "png", 300, 7)
myWidth <- myHeight*1.2

plotCex <- 1.2

vdHeight <- 7
vdWidth <- 7

signifThresh <- 0.05

source( file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg"), "set_dataset_colors.R"))
head(score_DT)
dataset_proc_colors <- setNames(score_DT$proc_col, score_DT$dataset)
length(dataset_proc_colors)

registerDoMC(ifelse(SSHFS, 2, 30))

# if(SSHFS) setwd("/media/electron/mnt/ed4/marie/scripts/Cancer_HiC_data_TAD_DA")
if(SSHFS) setwd("~/media/electron/mnt/ed4/marie/scripts/Cancer_HiC_data_TAD_DA")

outFold <- file.path("DATASETS_GENETAD_PVALS")
system(paste0("mkdir -p ", outFold))

logFile <- file.path(outFold, "dataset_genetad_pvals_logFile.txt")
system(paste0("rm -f ", logFile))

dsFold <- "PIPELINE/OUTPUT_FOLDER"
stopifnot(dir.exists(dsFold))

all_ds_files <- list.files(dsFold, recursive = TRUE, pattern="emp_pval_combined.Rdata", full.names = TRUE)

all_ds <- dirname(dirname(all_ds_files)) 
all_ds <- gsub(paste0(dsFold,"/"), "", all_ds)
stopifnot(length(all_ds) > 0)

txt <- paste0("... found # datasets:\t", length(all_ds), "\n")
printAndLog(txt, logFile)

txt <- paste0("... buildTable = ", buildTable, "\n")
printAndLog(txt, logFile)

txt <- paste0("... computeAUC = ", computeAUC, "\n")
printAndLog(txt, logFile)


txt <- paste0("... signifThresh = ", signifThresh, "\n")
printAndLog(txt, logFile)

##########################################################################################
##########################################################################################


curr_ds="ENCSR079VIJ_G401_40kb/TCGAkich_norm_kich"

# all_ds=all_ds[1:5]

if(buildTable){
  all_ds_geneTAD_pvals_list <- foreach(curr_ds = all_ds) %do% {
    # all_ds_DT <- foreach(curr_ds = topDS, .combine='rbind') %dopar% {
    txt <- paste0("*** START:\t", curr_ds, "\n")
    printAndLog(txt, logFile)

    ### RETRIEVE TAD GENES AND PVALUES
    cat("... retrieve TAD pvals \n")
    step11_fold <- file.path(dsFold, curr_ds, "11_runEmpPvalCombined")
    stopifnot(dir.exists(step11_fold))
    tadpvalFile <- file.path(step11_fold, "emp_pval_combined.Rdata")
    stopifnot(file.exists(tadpvalFile))
    tad_pval <- eval(parse(text = load(tadpvalFile)))
    tad_pval <- p.adjust(tad_pval, method = "BH")
    tadpval_DT <- data.frame(region=names(tad_pval), TAD_adjPval=tad_pval, stringsAsFactors = FALSE)
    rownames(tadpval_DT) <- NULL
    
    cat("... retrieve gene list \n")
    step0_fold <- file.path(dsFold, curr_ds, "0_prepGeneData")
    stopifnot(dir.exists(step0_fold))
    geneFile <- file.path(step0_fold, "pipeline_geneList.Rdata")
    stopifnot(file.exists(geneFile))
    geneList <- eval(parse(text = load(geneFile)))

    g2tFile <- file.path(dirname(curr_ds), "genes2tad", "all_genes_positions.txt")
    stopifnot(file.exists(g2tFile))
    g2t_dt <- read.delim(g2tFile, header=F, stringsAsFactors = FALSE, col.names=c( "entrezID", "chromo","start", "end", "region"))
    g2t_dt$entrezID <- as.character(g2t_dt$entrezID)
    
    g2t_dt <- g2t_dt[g2t_dt$entrezID %in% geneList,]
    stopifnot(nrow(g2t_dt) > 0)
    
    stopifnot(names(tad_pval) %in% g2t_dt$region)
    stopifnot(g2t_dt$region %in% names(tad_pval))
    
    g2t_dt$region <- as.character(g2t_dt$region)
    tadpval_DT$region <- as.character(tadpval_DT$region)
    g2tPval_DT <- merge(g2t_dt, tadpval_DT, by="region")
    
    gene_tad_pval <- setNames(g2tPval_DT$TAD_adjPval, as.character(g2tPval_DT$entrezID))
    
    gene_tad_pval

  } # end iterating datasets

  names(all_ds_geneTAD_pvals_list) <- all_ds
      
  outFile <-     file.path(outFold, "all_ds_geneTAD_pvals_list.Rdata")
  save(all_ds_geneTAD_pvals_list, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))
  
} else {# end if buildtable
    
  outFile <-     file.path(outFold, "all_ds_geneTAD_pvals_list.Rdata")
  all_ds_geneTAD_pvals_list <- eval(parse(text = load(outFile)))
  
}
# stop("-- ok \n")
# load("DATASETS_GENETAD_PVALS/all_ds_geneTAD_pvals_list.Rdata")

commonGenes <- Reduce(intersect, lapply(all_ds_geneTAD_pvals_list, names))

all_ds_geneTAD_pvals_commonList <- lapply(all_ds_geneTAD_pvals_list, function(x) {
  stopifnot(commonGenes %in% names(x) )
  x[commonGenes]
})

all_ds_geneTAD_pvals_DT <- as.data.frame(all_ds_geneTAD_pvals_commonList)
all_ds_geneTAD_pvals_DT <- t(all_ds_geneTAD_pvals_DT)
stopifnot(!any(is.na(all_ds_geneTAD_pvals_DT)))

# for each of the TADs -> in how many dataset it is signif
nSignif_byGenes <- apply(all_ds_geneTAD_pvals_DT, 2, function(x) sum(na.omit(x) <= signifThresh ))
range(nSignif_byGenes)

myTit <- paste0("# datasets in which TAD signif.")
myxlab <- "# datasets in which signif."
myylab <- "density"
mySub <- paste0("(signif. threshold = ", signifThresh, ")")
top5 <- names(sort(nSignif_byGenes, decreasing=T))[1:5]

legTxt <-  paste0("# datasets = ", nrow(all_ds_geneTAD_pvals_DT),"\n", "# TADs = ", ncol(all_ds_geneTAD_pvals_DT), "\n", "top5:\n", paste0(top5, collapse=",\n"))

outFile <- file.path(outFold, paste0("nbr_all_ds_in_which_signif.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot(density(nSignif_byGenes),
     main=myTit,
     xlab=myxlab,
     ylab=myylab,
     cex.axis=plotCex, cex.lab=plotCex
     )
mtext(side=3, text = mySub)
legend( "topright", bty="n",legend=legTxt, lty=-1)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


rankMax <- nrow(all_ds_geneTAD_pvals_DT)

rankVect <- seq(from=0, to=1, by=0.01)

cat(paste0("... start computing auc\n"))
if(computeAUC) {
  all_TADs_auc <- foreach(i = 1:ncol(all_ds_geneTAD_pvals_DT), .combine='c') %dopar% {
    curr_genetad_pvals <- all_ds_geneTAD_pvals_DT[,i]
    # rankVect <- 1:rankMax
    cumsum_ranks <- sapply(rankVect, function(x) sum(na.omit(curr_genetad_pvals) <= x  ))
    curr_auc <- auc(x = rankVect, y = cumsum_ranks)
    curr_auc
  }
  names(all_TADs_auc) <- colnames(all_ds_geneTAD_pvals_DT)
  outFile <- file.path(outFold, "all_geneTADs_auc.Rdata")
  save(all_TADs_auc, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))
} else {
  outFile <- file.path(outFold, "all_geneTADs_auc.Rdata")
  all_TADs_auc <- eval(parse(text = load(outFile)))
}

# 

entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
entrezDT <- read.delim(entrezDT_file, header=TRUE, stringsAsFactors = FALSE)
entrezDT$entrezID <- as.character(entrezDT$entrezID)

nTopGenes=5
topGenes <- sort(all_TADs_auc, decreasing=TRUE)[1:nTopGenes]
topGenesNames <- as.character(unlist(sapply(names(topGenes), function(x) entrezDT$symbol[as.character(entrezDT$entrezID) == x])))



cat(paste0("... start drawing\n"))
outFile <- file.path(outFold, paste0("all_geneTADs_cumsum_linePlot.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot(NULL,
     xlim=c(0,1),
     ylim=c(1,nrow(all_ds_geneTAD_pvals_DT)),
     main=paste0("cumsum # datasets <= pval thresh. (",signifThresh , ")"),
     ylab = paste0("# datasets"),
     xlab = paste0("gene TAD pvals")
     )
for(i in 1:ncol(all_ds_geneTAD_pvals_DT)) {
  curr_genetad_pvals <- all_ds_geneTAD_pvals_DT[,i]
  # rankVect <- 0:1
  cumsum_ranks <- sapply(rankVect, function(x) sum(na.omit(curr_genetad_pvals) <= x  ))
  lines(x=rankVect, y= cumsum_ranks)
}

legend("bottomright", topGenesNames, bty="n")

foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))



######################################################################################
######################################################################################
######################################################################################
cat(paste0("... written: ", logFile, "\n"))
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))

