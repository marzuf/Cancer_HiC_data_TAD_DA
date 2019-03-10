startTime <- Sys.time()
cat(paste0("> Rscript scores_vs_other_variables_withBoxplots.R\n"))

#  Rscript scores_vs_other_variables_withBoxplots.R

library(foreach)
library(doMC)

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")
setDir <- ifelse(SSHFS, "~/media/electron", "")

buildTable <- TRUE

rangeOffset <- 0.15
source(file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_11_18/analysis_utils.R")))
source(file.path(setDir, "/mnt/ed4/marie/scripts/EZH2_final_MAPQ/ezh2_utils_fct.R"))
source( file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg"), "set_dataset_colors.R"))
head(score_DT)
source( file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_11_18"), "analysis_utils.R"))
source( file.path("colors_utils.R"))
dataset_proc_colors <- setNames(score_DT$proc_col, score_DT$dataset)
length(dataset_proc_colors)

require(ggplot2)
source( file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg"), "set_dataset_colors.R"))
head(score_DT)
source( file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_11_18"), "analysis_utils.R"))
source( file.path("colors_utils.R"))
dataset_proc_colors <- setNames(score_DT$proc_col, score_DT$dataset)
length(dataset_proc_colors)

source("utils_fct.R")
source("utils_plot_fcts.R")
source("plot_lolliTAD_funct.R")

signifThresh <- 0.05
signifVar <- "adj.P.Val"

options(scipen=100)

registerDoMC(ifelse(SSHFS, 2, 40))

printAndLog <- function(txt, logfile) {
  cat(txt)
  cat(txt, file = logfile, append=T)
}

script0_name <- "0_prepGeneData"
script1_name <- "1_runGeneDE"
script11_name <- "11_runEmpPvalCombined"
script17_name <- "170revision2EZH2_score_auc_pval_permGenes"

famType1 <- "hgnc"
famType2 <- "hgnc_family_short"

plotType <- "svg"
myHeight <- ifelse(plotType == "png", 400, 10)
# myWidth <- ifelse(plotType == "png", 600, 10)
myWidth <- myHeight


plotAxisCex <- 1.2

aucCoexprDistFolder <- file.path("AUC_COEXPRDIST_WITHFAM_SORTNODUP")
stopifnot(dir.exists(aucCoexprDistFolder))

pipOutFolder <- file.path("PIPELINE", "OUTPUT_FOLDER")

stopifnot(file.exists(pipOutFolder))
settingFilesFolder <- file.path("PIPELINE", "INPUT_FILES")

stopifnot(file.exists(settingFilesFolder))

geneVarFile <- "GENE_VARIANCE/LOG2FPKM/all_ds_geneVarDT.Rdata"

stopifnot(file.exists(geneVarFile))
cat(paste0("!!! variance retrieved from: ", geneVarFile, "\n"))

all_datasets <- list.files(list.files(pipOutFolder, full.names = TRUE), full.names=TRUE)

all_datasets <- file.path(basename(dirname(all_datasets)), basename(all_datasets))

# hicds                exprds                                         dataset coexprDistAUC   fccAUC nSamp1
# 5 ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_wt_mutCTNNB1 ENCSR312KHQ_SK-MEL-5_40kb_TCGAskcm_wt_mutCTNNB1      1.151576 1.456819    351
# nSamp2 nAllSamp ratioSamp nGenes nSignifGenes nTADs nSignifTADs meanGenesByTAD medianGenesByTAD meanMostVar medianMostVar
# 5     12      363     29.25  13198          366  2368        2079        5.57348                5    3.147607      2.597039
# > 
stopifnot("ENCSR312KHQ_SK-MEL-5_40kb/TCGAskcm_wt_mutCTNNB1" %in% all_datasets)
# all_datasets <- all_datasets[basename(all_datasets) != "TCGAskcm_wt_mutCTNNB1"]

if("TCGAskcm_wt_mutCTNNB1" %in% basename(all_datasets)){ 
  outFold <- file.path("SCORES_VS_OTHER_VARIABLES_WITH_BOXPLOTS")
}else {
  outFold <- file.path("SCORES_VS_OTHER_VARIABLES_WITH_BOXPLOTS_NOOUT")
}
dir.create(outFold, recursive = TRUE)

logFile <- file.path(outFold, "score_vs_other_variables_logFile.txt")
if(!SSHFS) file.remove(logFile)
if(SSHFS) logFile <- ""

  
# head(all_datasets)
stopifnot(length(all_datasets) >  0 )
cat(paste0("... found ", length(all_datasets), " datasets\n"))

source( file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg"), "set_dataset_colors.R"))
head(score_DT)

dataset_proc_colors <- setNames(score_DT$proc_col, score_DT$dataset)
length(dataset_proc_colors)


curr_dataset = "NCI-H460_40kb/TCGAluad_nonsmoker_smoker"

# return NULL if file not found ??? [if yes -> build table skipping missing files, otherwise raise error and stop]
returnNull <-  FALSE

txt <- paste0("... geneVarFile\t=\t", geneVarFile, "\n")
printAndLog(txt, logFile)
txt <- paste0("... returnNull\t=\t", as.character(returnNull), "\n")
printAndLog(txt, logFile)
txt <- paste0("... found # ds\t=\t", length(all_datasets) , "\n" )
printAndLog(txt, logFile)

if(buildTable) {
  
  geneVarDT <- eval(parse(text = load(geneVarFile)))
  
  datasets_variables_DT <- foreach(curr_dataset = all_datasets, .combine='rbind') %dopar% {
    
    exprds <- basename(curr_dataset)
    hicds <- dirname(curr_dataset)
    
    curr_dataset_name <- paste0(hicds, "_", exprds)
    curr_dataset_path <- file.path(hicds, exprds)
    
    gene2tadDT_file <- file.path(hicds, "genes2tad", "all_genes_positions.txt")
    stopifnot(file.exists(gene2tadDT_file))
    gene2tad_DT <- read.delim(gene2tadDT_file, header=F, stringsAsFactors = F, col.names=c("entrezID", "chromo", "start", "end", "region"))
    gene2tad_DT$entrezID <- as.character(gene2tad_DT$entrezID)
    
    
    txt <- paste0("> START ", curr_dataset, "\n")
    printAndLog(txt, logFile)
    
    cat("... load auc coexpr dist\n")
    aucCoexprDistFile <- file.path(aucCoexprDistFolder, hicds, paste0(exprds, "_", famType1), famType2, "auc_values.Rdata")
    stopifnot(file.exists(aucCoexprDistFile))
    
    all_coexprDistAUC <- eval(parse(text = load(aucCoexprDistFile)))
    coexprDistAUC <- all_coexprDistAUC[["auc_ratio_same_over_diff_distVect"]]
    cat("coexprDistAUC = ", coexprDistAUC, "\n")
    stopifnot(is.numeric(coexprDistAUC))
    
    cat("... load auc FCC \n")
    step17_fold <- file.path(pipOutFolder, curr_dataset_path, script17_name)
    stopifnot(dir.exists(step17_fold))
    aucFCC_file <- file.path(step17_fold, "auc_ratios.Rdata")
    stopifnot(file.exists(aucFCC_file))
    

    all_fccAUC <- eval(parse(text = load(aucFCC_file)))
    
    fccAUC <- all_fccAUC["prodSignedRatio_auc_permGenes"]
    stopifnot(is.numeric(fccAUC))
    
    cat("... load pipeline_geneList\n")
    geneFile <- file.path(pipOutFolder, 
                          curr_dataset_path, 
                          script0_name,
                          "pipeline_geneList.Rdata")
    pipeline_geneList <- eval(parse(text = load(geneFile)))
    
    tmp_g2t <- gene2tad_DT[gene2tad_DT$entrezID %in% pipeline_geneList,]
    stopifnot(nrow(tmp_g2t) > 0)
    nGenesByTAD_dt <- data.frame(region = as.character(names(table(tmp_g2t$region))),
                                 nGenes = as.numeric(table(tmp_g2t$region)),
                                 stringsAsFactors = FALSE
                                 )
    meanGenesByTAD <- mean(nGenesByTAD_dt$nGenes)
    medianGenesByTAD <- median(nGenesByTAD_dt$nGenes)
    
    
    cat("... load limma DT\n")
    limmaFile <- file.path(pipOutFolder, 
                           curr_dataset_path, 
                           script1_name,
                           "DE_topTable.Rdata")

    limmaDT <- eval(parse(text = load(limmaFile)))
    
    stopifnot(limmaDT$genes == rownames(limmaDT))
    stopifnot(names(pipeline_geneList) %in% limmaDT$genes)
    
    limmaDT <- limmaDT[limmaDT$genes %in% names(pipeline_geneList),]
    stopifnot(nrow(limmaDT) > 0)

    cat("... load TAD\n")
    tadpvalFile <-  file.path(pipOutFolder, curr_dataset_path, script11_name, "emp_pval_combined.Rdata")


    tad_pval <- eval(parse(text = load(tadpvalFile)))
    tad_pval <- p.adjust(tad_pval, method = "BH")
    
    
    cat("... load setting files\n")
    # PIPELINE/INPUT_FILES/NCI-H460_40kb/run_settings_TCGAluad_mutKRAS_mutEGFR.R
    
    settingF <- file.path(settingFilesFolder,
                          hicds,
                          paste0("run_settings_", exprds, ".R"))

    source(settingF)
    
    sample1_file <- file.path(setDir, sample1_file)
    s1 <- eval(parse(text = load(sample1_file)))
    
    sample2_file <- file.path(setDir, sample2_file)
    s2 <- eval(parse(text = load(sample2_file)))
    
    nGenes <- length(pipeline_geneList)
    nSignifGenes <- sum(limmaDT[,signifVar] <= signifThresh)
    
    nTADs <- length(tad_pval)
    nSignifTADs <- sum(tad_pval <= signifThresh)
    
    nSamp1 <- length(s1)
    nSamp2 <- length(s2)
    
    stopifnot(curr_dataset_name %in% geneVarDT$dataset)
    meanMostVar <- geneVarDT$meanMostVar[geneVarDT$dataset == curr_dataset_name]
    medianMostVar <- geneVarDT$medianMostVar[geneVarDT$dataset == curr_dataset_name]
    

    data.frame(
      hicds = hicds,
      exprds = exprds,
      dataset=curr_dataset_name,
      coexprDistAUC = coexprDistAUC,
      fccAUC = fccAUC,
      nSamp1 = nSamp1,
      nSamp2 = nSamp2,
      nAllSamp = nSamp1+nSamp2,
      ratioSamp = nSamp1/nSamp2,
      log2ratioSamp = log2(nSamp1/nSamp2),
      nGenes = nGenes,
      nSignifGenes = nSignifGenes,
      nTADs = nTADs,
      nSignifTADs=nSignifTADs,
      meanGenesByTAD = meanGenesByTAD,
      medianGenesByTAD = medianGenesByTAD,
      meanMostVar = meanMostVar,
      medianMostVar = medianMostVar,
      stringsAsFactors = FALSE
    )
  }
  rownames(datasets_variables_DT) <- NULL
  stopifnot(nrow(datasets_variables_DT) == length(all_datasets))
  outFile <- file.path(outFold, "datasets_variables_DT.Rdata")
  save(datasets_variables_DT, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))
  
} else {
  outFile <- file.path(outFold, "datasets_variables_DT.Rdata")
  #outFile <- "SCORES_VS_OTHER_VARIABLES_WITH_BOXPLOTS/datasets_variables_DT.Rdata"
  stopifnot(file.exists(outFile))
  datasets_variables_DT <- eval(parse(text = load(outFile)))
}

noVar <- which(apply(datasets_variables_DT, 2, function(x) length(unique(x) == 1)) == 1)

if(length(noVar) > 0) {
  cat("... unique value for: ", paste0(colnames(datasets_variables_DT)[noVar], collapse=";"), "\n")
  datasets_variables_DT[,noVar] <- NULL
}
# => remove median -> same single value -> cannot add the curve to the plot

# load("SCORES_VS_OTHER_VARIABLES_WITH_BOXPLOTS/datasets_variables_DT.Rdata")

# put it before adding color and type columns...
all_vars <- colnames(datasets_variables_DT)[! colnames(datasets_variables_DT) %in% c("exprds", "hicds", "dataset", "nSamp1", "nSamp2")]
all_vars <- colnames(datasets_variables_DT)[! colnames(datasets_variables_DT) %in% c("exprds", "hicds", "dataset")]


datasets_variables_DT$color <- unlist(sapply(datasets_variables_DT$exprds, function(x) cancer_subColors[cancer_subAnnot[x]]))
stopifnot(!is.na(datasets_variables_DT$color))

datasets_variables_DT$dsType <- unlist(sapply(datasets_variables_DT$exprds, function(x) cancer_subAnnot[x]))
stopifnot(!is.na(datasets_variables_DT$dsType))

titNames <- c(
  #  coexprDistAUC="AUC ratio - pairwise coexpr.",
  #  fccAUC = "AUC ratio - FCC",
  coexprDistAUC="% AUC ratio increase - pairwise coexpr.",
  fccAUC = "% AUC ratio increase - FCC",
  nSamp1 = "# samples (cond1)",
  nSamp2 = "# samples (cond2)",
  log2ratioSamp = "lgo2ratio cond1/cond2 samples",
  nAllSamp = "# samples (all)",
  ratioSamp = "ratio cond1/cond2 samples",
  nGenes = "# genes",
  nSignifGenes = paste0("# signif. DE genes (", signifVar, " <= ", signifThresh, ")"),
  nTADs = "# TADs",
  nSignifTADs = paste0("# signif. DA TADs (", signifVar, " <= ", signifThresh, ")"),
  meanGenesByTAD = "mean # genes/TAD",
  medianGenesByTAD = "median # genes/TAD",
  medianMostVar = "median var. of most var. genes",
  meanMostVar = "mean var. of most var. genes"
)

######################################################################################
###################################################################################### # draw the boxplots and multidens
######################################################################################

curr_var  = "nAllSamp"

datasets_variables_DT$dsTypeLabel <- datasets_variables_DT$dsType
stopifnot(!is.na(datasets_variables_DT$dsTypeLabel))

# all_vars="medianGenesByTAD"
head(datasets_variables_DT)

stopifnot(all_vars %in% names(titNames))

for(curr_var in all_vars) {
  outFile <- file.path(outFold, paste0(curr_var, "_", "cancer_TCGA_boxplot.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  boxplot(as.formula(paste0(curr_var, " ~  dsTypeLabel")), 
          data = datasets_variables_DT,
          boxfill = cancer_subColors[as.character(levels(as.factor(datasets_variables_DT$dsTypeLabel)))],
          main = paste0(curr_var),
#          ylab = paste0(curr_var),
          ylab = paste0(titNames[curr_var]),
          las = 2,
           cex.lab = plotAxisCex, cex.axis = plotAxisCex
  )  
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  outFile <- file.path(outFold, paste0(curr_var, "_", "cancer_TCGA_multidens.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  plot_multiDens_setcols(
    lapply(split(datasets_variables_DT, datasets_variables_DT$dsType), function(x) {
      tmp <- x[[curr_var]]
      # if(mytransf == "log10") tmp <- log10(tmp)
      tmp
    }),
    my_cols = cancer_subColors[as.character(levels(as.factor(datasets_variables_DT$dsType)))],
    my_xlab = paste0(titNames[curr_var]),
    plotTit =  paste0(curr_var)
  )
  # mtext(side=3, text=subTit)
  foo <- dev.off() 
  cat(paste0("... written: ", outFile, "\n"))
}

######################################################################################
######################################################################################

stopifnot(!duplicated(datasets_variables_DT$dataset))
rownames(datasets_variables_DT) <- datasets_variables_DT$dataset
datasets_variables_DT$dataset <- NULL

other_vars <- all_vars

offSets <- c(
#  coexprDistAUC=0.03,
#  fccAUC = 0.03,
  coexprDistAUC=0.03*100,
  fccAUC = 0.03*100,
  nSamp1 = 150,
  nSamp2 = 150,
  nAllSamp = 150,
  ratioSamp = 5,
  log2ratioSamp = log2(5),
  nGenes = 1000,
  nSignifGenes = 1000,
  nTADs = 100,
  nSignifTADs = 500,
  meanGenesByTAD = 2,
  medianGenesByTAD = 2,
  meanMostVar = 0.1,
  medianMostVar = 0.1
)

# stopifnot(colnames(datasets_variables_DT) %in% names(titNames))

mySub <- paste0("all datasets (n=", nrow(datasets_variables_DT), ")")

curr_colors <- as.character(cancer_subColors[as.character(cancer_subAnnot[as.character(datasets_variables_DT$exprds)])])

stopifnot(!is.na(curr_colors))
my_colors <- cancer_subColors

for(ref_var in c("coexprDistAUC", "fccAUC")) {
  for(curr_var in other_vars ) {
    
    outFile <- file.path(outFold, paste0(ref_var, "_", curr_var, ".", plotType))
    
    myx <- datasets_variables_DT[,curr_var]
    myy <- datasets_variables_DT[,ref_var]

    if(curr_var == "coexprDistAUC" | curr_var == "fccAUC") myx <- (myx-1) * 100
    if(ref_var == "coexprDistAUC" | ref_var == "fccAUC") myy <- (myy-1) * 100
    
    do.call(plotType, list(outFile, height = myHeight, width = myWidth))
    my_plot_function(varX = myx, varY=myy, 
                     mylabels = as.character(rownames(datasets_variables_DT)), 
                     withLab=F,
                     main = paste0(titNames[ref_var], " vs. ", titNames[curr_var]),
                     xlab = paste0(titNames[curr_var]),
                     ylab = paste0(titNames[ref_var]),
                     xlim = range(myx) + c(-offSets[curr_var], offSets[curr_var]),
                     ylim = range(myy) + c(-offSets[ref_var], offSets[ref_var]),
                     cex.lab = plotAxisCex, cex.axis = plotAxisCex
                     )
    text(x = myx, y = myy,
         labels =  as.character(rownames(datasets_variables_DT)), 
         pch=16,
         col = curr_colors,
         bty="n",
         cex=0.7)
    mtext(side=3, text = mySub)
    legend("bottomright",
           legend=names(my_colors),
           lty=1,
           col = my_colors,
           lwd = 5,
           bty="n",
           cex = 0.7)
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
    if(grepl("^n.*Genes", curr_var) | grepl("TADs", curr_var)) {
      myx <- datasets_variables_DT[,curr_var]
      myy <- datasets_variables_DT[,ref_var]
      
      myx <- log10(myx+1)
      
      if(grepl("^n.*Genes", curr_var)) {
      xoffset <- ifelse(curr_var == "nGenes", 0.05,
                        ifelse(curr_var == "nSignifGenes", 0.35, stop("")))
      } else if(grepl("TADs", curr_var)) {
        xoffset <- ifelse(curr_var == "nTADs", 0.05,
                          ifelse(curr_var == "nSignifTADs", 0.35, stop("")))
      } else {stop("")}
      outFile <- file.path(outFold, paste0(ref_var, "_", curr_var, "_log10.", plotType))
      
      do.call(plotType, list(outFile, height = myHeight, width = myWidth))
      my_plot_function(varX = myx, varY=myy, 
                       mylabels = as.character(rownames(datasets_variables_DT)), 
                       withLab=F,
                       main = paste0(titNames[ref_var], " vs. ", titNames[curr_var]),
                       xlab = paste0(titNames[curr_var], "(log10[#+1])" ),
                       ylab = paste0(titNames[ref_var]),
                       xlim = range(myx) + c(-xoffset, xoffset),
                       ylim = range(myy) + c(-offSets[ref_var], offSets[ref_var]),
                       cex.lab = plotAxisCex, cex.axis = plotAxisCex
      )
      text(x = myx, y = myy,
           labels =  as.character(rownames(datasets_variables_DT)), 
           pch=16,
           col = curr_colors,
           bty="n",
           cex=0.7)
      mtext(side=3, text = mySub)
      legend("bottomright",
             legend=names(my_colors),
             lty=1,
             col = my_colors,
             lwd = 5,
             bty="n",
             cex = 0.7)
      foo <- dev.off()
      cat(paste0("... written: ", outFile, "\n"))
      
    }
  } # end iterating over variables
} # end iterating FCC - coexprdist

all_cmps <- list(
c(curr_var = "nSignifGenes", ref_var= "meanMostVar"),
c(curr_var = "nSignifTADs", ref_var= "meanMostVar"),

c(curr_var = "nSignifGenes", ref_var= "medianMostVar"),
c(curr_var = "nSignifTADs", ref_var= "medianMostVar"),

c(curr_var = "nSignifGenes", ref_var= "nSignifTADs")
)

for(i in seq_along(all_cmps)){
    curr_var = all_cmps[[i]]["curr_var"]
    ref_var = all_cmps[[i]]["ref_var"]

    myx <- log10(datasets_variables_DT[, curr_var]+1)
    myy <- datasets_variables_DT[, ref_var]

    outFile <- file.path(outFold, paste0(ref_var, "_", curr_var, "_log10.", plotType))
    do.call(plotType, list(outFile, height = myHeight, width = myWidth))
    my_plot_function(varX = myx, varY=myy, 
                     mylabels = as.character(rownames(datasets_variables_DT)), 
                     withLab=F,
                     main = paste0(titNames[ref_var], " vs. ", titNames[curr_var]),
                     xlab = paste0(titNames[curr_var], "(log10[#+1])" ),
                     ylab = paste0(titNames[ref_var]),
                     xlim = range(myx) + c(-1,1),# + c(-log10(xoffset[curr_var]), log10(xoffset[curr_var])),
                     ylim = range(myy) + c(-offSets[ref_var], offSets[ref_var]),
                     cex.lab = plotAxisCex, cex.axis = plotAxisCex
    )
    text(x = myx, y = myy,
         labels =  as.character(rownames(datasets_variables_DT)), 
         pch=16,
         col = curr_colors,
         bty="n",
         cex=0.7)
    mtext(side=3, text = mySub)
    legend("bottomright",
           legend=names(my_colors),
           lty=1,
           col = my_colors,
           lwd = 5,
           bty="n",
           cex = 0.7)
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))

    myx <- datasets_variables_DT[, curr_var]
    myy <- datasets_variables_DT[, ref_var]

    outFile <- file.path(outFold, paste0(ref_var, "_", curr_var, ".", plotType))
    do.call(plotType, list(outFile, height = myHeight, width = myWidth))
    my_plot_function(varX = myx, varY=myy, 
                     mylabels = as.character(rownames(datasets_variables_DT)), 
                     withLab=F,
                     main = paste0(titNames[ref_var], " vs. ", titNames[curr_var]),
                     xlab = paste0(titNames[curr_var], "" ),
                     ylab = paste0(titNames[ref_var]),
                     xlim = range(myx) + c(-1000,1000),#+ c(-xoffset[curr_var], xoffset[curr_var]),
                     ylim = range(myy) + c(-offSets[ref_var], offSets[ref_var]),
                     cex.lab = plotAxisCex, cex.axis = plotAxisCex
    )
    text(x = myx, y = myy,
         labels =  as.character(rownames(datasets_variables_DT)), 
         pch=16,
         col = curr_colors,
         bty="n",
         cex=0.7)
    mtext(side=3, text = mySub)
    legend("bottomright",
           legend=names(my_colors),
           lty=1,
           col = my_colors,
           lwd = 5,
           bty="n",
           cex = 0.7)
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))

}

######################################################################################
######################################################################################
######################################################################################
cat(paste0("... written: ", logFile, "\n"))
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))
