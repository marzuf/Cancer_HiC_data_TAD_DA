
# Rscript pval_by_subtypes.R # default 0.05
# Rscript pval_by_subtypes.R 0.05

startTime <- Sys.time()

cat("> START pval_by_subtypes.R \n")

SSHFS <- FALSE

require(foreach)
require(doMC)

source("utils_fct.R")
source("utils_plot_fcts.R")

setDir <- ifelse(SSHFS, "~/media/electron", "")


source( file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg"), "set_dataset_colors.R"))
head(score_DT)
source( file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_11_18"), "analysis_utils.R"))
source( file.path("colors_utils.R"))
dataset_proc_colors <- setNames(score_DT$proc_col, score_DT$dataset)
length(dataset_proc_colors)


# head(cancer_subAnnot)
# TCGAstad_msi_gs            GSE79209_dysp_nodysp GSE58135_tripleNeg_adjTripleNeg 
# "subtypes"                       "lesions"                     "vs_normal" 
# head(dataset_proc_colors)
# TCGAgbm_classical_proneural              TCGAgbm_classical_mesenchymal 
# "green4"                                   "green4" 
# TCGAskcm_lowInf_highInf                           TCGAcoad_msi_mss 
# "green4"                                   "green4"

plotType <- "svg"
myHeight <- ifelse(plotType=="png", 500, 7)
myWidth <- ifelse(plotType=="png", 600, 9)


registerDoMC(ifelse(SSHFS, 2, 40))

build_allPvals_list <- F


signifThresh <- 0.05

args <- commandArgs(trailingOnly = TRUE)
if(length(args) > 0) {
  signifThresh <- as.numeric(args[1])
  stopifnot(!is.na(signifThresh))
}

stopifnot(signifThresh >= 0)
stopifnot(signifThresh <= 1)

outFolder <- file.path("PVAL_BY_SUBTYPES", paste0("signif", signifThresh))
dir.create(outFolder, recursive = TRUE)

logFile=""

stopifnot(!is.na(signifThresh))
stopifnot(signifThresh > 0)

#PIPELINE/OUTPUT_FOLDER/GSE105318_DLD1_40kb/TCGAcoad_msi_mss//emp_pval_combined.Rdata
script0_name <- "0_prepGeneData"
script1_name <- "1_runGeneDE"
script9_name <- "9_runEmpPvalMeanTADLogFC"
script10_name <- "10_runEmpPvalMeanTADCorr"
script11_name <- "11_runEmpPvalCombined"

pipOutFolder <- file.path("PIPELINE", "OUTPUT_FOLDER")

all_hicexpr_ds <- unname(unlist(sapply(list.files(pipOutFolder, full.names = TRUE), function(x) file.path(basename(x),  list.files(x)))))
stopifnot(dir.exists(file.path(pipOutFolder, all_hicexpr_ds)))

all_exprds <- data.frame(exprds=basename(all_hicexpr_ds), stringsAsFactors = FALSE)
subTypeDT <- data.frame(exprds=names(cancer_subAnnot), subtype=cancer_subAnnot, stringsAsFactors = FALSE)
colorDT <- data.frame(exprds=names(dataset_proc_colors), color=dataset_proc_colors, stringsAsFactors = FALSE)
stopifnot(all_exprds$exprds %in% subTypeDT$exprds)
stopifnot(all_exprds$exprds %in% colorDT$exprds)

subtype_data_DT <- unique(merge(merge(all_exprds, subTypeDT, by="exprds", all.x=T, all.y=F), 
                         colorDT, by="exprds", all.x=T, all.y=F))
stopifnot(!is.na(subtype_data_DT))

ds=all_hicexpr_ds[1]

if(build_allPvals_list){
  
  cat("... start building allPvals_allDS_DT data \n")
  
  allPvals_list <- foreach(ds = all_hicexpr_ds) %dopar% {
    
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
    
    # RETRIEVE THE GENES PVALS
    stopifnot(dir.exists(file.path(dsPipOutDir, script0_name)))
    gene_listFile <- file.path(dsPipOutDir, script0_name, "pipeline_geneList.Rdata")
    stopifnot(file.exists(gene_listFile))
    geneList <- eval(parse(text = load(gene_listFile)))
    # ENSG00000000003 ENSG00000000005 ENSG00000000419 ENSG00000000457 ENSG00000000460 
    # "7105"         "64102"          "8813"         "57147"         "55732" 
    
    stopifnot(dir.exists(file.path(dsPipOutDir, script1_name)))
    gene_pvalFile <- file.path(dsPipOutDir, script1_name, "DE_topTable.Rdata")
    stopifnot(file.exists(gene_pvalFile))
    topTable_DT <- eval(parse(text = load(gene_pvalFile)))
    # genes      logFC     AveExpr         t      P.Value
    # ENSG00000116745 ENSG00000116745 -1.2849837  0.05196191 -3.895238 0.0002779426
    
    stopifnot(names(geneList) %in% topTable_DT[, "genes"])
    gene_pval <- setNames(topTable_DT[, "adj.P.Val"], topTable_DT[, "genes"])
    
    gene_pval <- gene_pval[names(geneList)]
    tmp <- names(gene_pval)
    stopifnot(length(gene_pval) == length(geneList))
    names(gene_pval) <- as.character(geneList)
    stopifnot(tmp == names(geneList))
    
    gene_pval <- sort(gene_pval)
    
    list(
      hicds=hicds,
      exprds=exprds,
      tad=all_tads,
      tad_id=paste(hicds, exprds, all_tads, sep="_"),
      adj_pvalFC = tad_pvalFC[all_tads],
      adj_pvalCorr = tad_pvalCorr[all_tads],
      adj_pvalComb = tad_pvalComb[all_tads],
      gene_pvals = gene_pval,
      stringsAsFactors = FALSE
    )

  }
  names(allPvals_list) <- all_hicexpr_ds
  outFile <- file.path(outFolder, "allPvals_list.Rdata")
  save(allPvals_list, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))
  
} else { # if not to build
  outFile <- file.path(outFolder, "allPvals_list.Rdata")
  allPvals_list <- eval(parse(text = load(outFile)))
}


# load("PVAL_BY_SUBTYPES/signif0.05/allPvals_list.Rdata")

# head(allPvals_list)
# str(allPvals_list)

# > names(allPvals_list[[1]])
# [1] "hicds"            "exprds"           "tad"              "tad_id"           "adj_pvalFC"      
# [6] "adj_pvalCorr"     "adj_pvalComb"     "gene_pvals"       "stringsAsFactors"

vars_densplot <- c("adj_pvalFC", "adj_pvalCorr", "adj_pvalComb", "gene_pvals")

my_subtypes_name <- unique(subtype_data_DT$subtype)
my_subtypes_color <- unique(subtype_data_DT$subtype)


var="adj_pvalFC"
foo <- foreach(var = vars_densplot) %dopar% {
  
  var_allds_DT <- do.call('rbind', lapply(seq_along(allPvals_list), function(x) {
   tmp <- data.frame(exprds=basename(names(allPvals_list)[x]), pval= allPvals_list[[x]][[var]], row.names = NULL)
   colnames(tmp) <- c("exprds", var)
   tmp
  } ))
  
  var_allds_DT_with_annot <- merge(var_allds_DT, subtype_data_DT, by="exprds", all.x=T, all.y=F)
  stopifnot(!is.na(var_allds_DT_with_annot))
  
  my_xlab <- paste0(var)
  
  dsByType <- table(unique(var_allds_DT_with_annot[, c("exprds", "subtype")])$subtype)
  nDSbyType <- setNames(as.numeric(dsByType), names(dsByType))
  
  subTit <- paste0("# DS: ", paste0(names(nDSbyType), "=", as.numeric(nDSbyType), collapse = " - "))
  
  cropTo <- 5
  
  for(tokeep in c("all", "signif")) {
    
    myTit <- paste0(var)
    if(tokeep == "signif") myTit <- paste0(myTit, " - ", tokeep, " (p-val <= ", signifThresh, ")")
    
    outFile <- file.path(outFolder, paste0("multidens_bySubtypes_", var , "_", tokeep, ".", plotType))
    do.call(plotType, list(outFile, height=myHeight, width=myWidth))
    plot_multiDens_setcols(
      lapply(split(var_allds_DT_with_annot, var_allds_DT_with_annot$subtype), function(x) {
        my_values <- x[[var]]
        my_values[my_values < 10^-cropTo] <- 10^-cropTo
        if(tokeep == "signif") my_values <- my_values[my_values <= signifThresh]
        my_values
      }),
      my_cols = cancer_subColors[as.character(levels(as.factor(var_allds_DT_with_annot$subtype)))],
      my_xlab = paste0(var),
      plotTit = myTit
    )
    mtext(side=3, text=subTit)
    foo <- dev.off()  
    cat(paste0("... written: ", outFile, "\n"))
    
    outFile <- file.path(outFolder, paste0("multidens_bySubtypes_", var , "_", tokeep, "_log10.", plotType))
    do.call(plotType, list(outFile, height=myHeight, width=myWidth))
    plot_multiDens_setcols(
      lapply(split(var_allds_DT_with_annot, var_allds_DT_with_annot$subtype), function(x) {
        my_values <- log10(x[[var]])
        my_values[my_values < -cropTo] <- (-cropTo)
        if(tokeep == "signif") my_values <- my_values[my_values <= log10(signifThresh)]
        my_values }),
      my_cols = cancer_subColors[as.character(levels(as.factor(var_allds_DT_with_annot$subtype)))],
      my_xlab = paste0(var, "[log10]"),
      legPos = "topleft",
      plotTit = myTit
    )
    mtext(side=3, text=subTit)
    foo <- dev.off()  
    cat(paste0("... written: ", outFile, "\n"))
    
    
  }
  
  
  
  
  
}


# xvar <- "adj_pvalFC"
# yvar <- "adj_pvalCorr"
# myTit <- "all TCGA datasets - all TADs"
# mySub <- paste0("(# TADs = ", length(allPvals_allDS_DT$tad_id) , ")")
# 
# plotCex <- 1.2

# outFile <- file.path(outFolder, paste0("pvalCorr_vs_pvalFC_allTADs", ".", plotType))
# do.call(plotType, list(outFile, height=myHeight, width=myWidth))
# myx <- allPvals_allDS_DT[,xvar]
# myy <- allPvals_allDS_DT[,yvar]
# densplot(x = myx,
#          xlab=paste0(xvar),
#          xlim=range(myx,na.rm=T),
#          ylim=range(myy,na.rm=T),
#          y = myy,
#          ylab=paste0(yvar),
#          main = myTit,
#          cex=0.7,
#          cex.axis=plotCex,
#          cex.lab=plotCex,
#          pch=16
# )
# mtext(side=3, text=mySub)
# foo <- dev.off()
# cat(paste0("... written: ", outFile, "\n"))

# plot_multiDens(list(
#   all_matchingRatio = all_matchDT$matchingRatio,
#   best_matchingRatio = all_bestMatchDT$matchingRatio),
#   plotTit="matchingRatio", legTxt=NULL, legPos="topleft", my_ylab="density", my_xlab=""
# )



# ######################################################################################
# ######################################################################################
# ######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))
