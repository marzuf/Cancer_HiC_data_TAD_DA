
# Rscript coexpr_curves.R

script_name <- "coexpr_curves.R"

startTime <- Sys.time()

cat("> START coexpr_curves.R \n")

SSHFS <- FALSE

buildData <- TRUE

require(foreach)
require(doMC)
registerDoMC(ifelse(SSHFS, 2, 40))

# source("utils_fct.R")

plotType <- "png"
myHeight <- ifelse(plotType=="png", 500, 7)
myWidth <- myHeight
plotCex <- 1.4


outFolder <- "COEXPR_CURVES"
dir.create(outFolder, recursive=TRUE)

corrMet <- "pearson"
script0_name <- "0_prepGeneData"
mypattern <- "coexprDT.Rdata"

coexprFiles  <- list.files("CREATE_COEXPR_SORTNODUP", recursive = TRUE, pattern = mypattern)
stopifnot(length(coexprFiles) > 0)


coexprDataFile = coexprFiles[1]
hicds = "K562_40kb"
exprds = "TCGAlaml_wt_mutFLT3"
coexprDataFile = file.path("CREATE_COEXPR_SORTNODUP", 
                            hicds, 
                            paste0(exprds, "_", corrMet),
                            "coexprDT.Rdata")



coexprFilesOutNames <- coexprFiles
coexprFilesOutNames <- gsub(file.path(paste0("_", corrMet), paste0(mypattern)),  "", coexprFilesOutNames)

if(buildData) {

  allData_within_between_coexpr <- foreach(coexprDataFile = coexprFiles) %do% {
    
    cat("... start ", coexprDataFile, " \n")
    
    
    hicds <- dirname(dirname(coexprDataFile))
    stopifnot(dir.exists(hicds))
    exprds <- basename(dirname(coexprDataFile))
    exprds <- gsub(paste0("_", corrMet), "", exprds)
    
    cat("... start ", hicds, " - ", exprds, "\n")
    
    coexprDataFile = file.path("CREATE_COEXPR_SORTNODUP", 
                               hicds, 
                               paste0(exprds, "_", corrMet),
                               "coexprDT.Rdata")
    
    
    cat(coexprDataFile, "\n")
    # cat(coexprFiles[1], "\n")
    stopifnot(file.exists(coexprDataFile))
    
    cat("... load coexprDT...\n")
    
    coexprDT <- eval(parse(text = load(coexprDataFile)))
    head(coexprDT)
    coexprDT$gene1 <- as.character(coexprDT$gene1)
    coexprDT$gene2 <- as.character(coexprDT$gene2)
    
    g2tFile <- file.path(hicds, "genes2tad", "all_genes_positions.txt")
    stopifnot(file.exists(g2tFile))
    g2t_DT <- read.delim(g2tFile, header=F, 
                         col.names = c("entrezID",  "chromo", "start", "end", "region"), stringsAsFactors = FALSE)
    g2t_DT$entrezID <- as.character(g2t_DT$entrezID)
    
    dsPipOutDir <- file.path("PIPELINE", "OUTPUT_FOLDER", hicds, exprds)
    stopifnot(dir.exists(dsPipOutDir))
    stopifnot(dir.exists(file.path(dsPipOutDir, script0_name)))
    geneListFile <- file.path(dsPipOutDir, script0_name, "pipeline_geneList.Rdata")
    stopifnot(file.exists(geneListFile))
    pipeline_geneList <- eval(parse(text = load(geneListFile))) # 
    
    regionListFile <- file.path(dsPipOutDir, script0_name, "pipeline_regionList.Rdata")
    stopifnot(file.exists(regionListFile))
    pipeline_regionList <- eval(parse(text = load(regionListFile))) # 
    
    stopifnot(pipeline_geneList %in% g2t_DT$entrezID)
    g2t_DT <- g2t_DT[g2t_DT$entrezID %in% pipeline_geneList & g2t_DT$region %in% pipeline_regionList, ]
    stopifnot(g2t_DT$entrezID %in% coexprDT$gene1 | g2t_DT$entrezID %in% coexprDT$gene2 )
    stopifnot(pipeline_regionList %in% g2t_DT$region)
    
    reg = pipeline_regionList[1]
    within_between_coexpr_data <- foreach(reg = pipeline_regionList) %dopar% {
      
      cat(paste0("...... start ", hicds, " - ", exprds, " - ", reg, " \n"))
      
      tad_genes <- g2t_DT$entrezID[g2t_DT$region == reg]
      tad_coexprDT <- coexprDT[coexprDT$gene1 %in% tad_genes | coexprDT$gene2 %in% tad_genes,]
    
      stopifnot(nrow(tad_coexprDT) > 0)  
      stopifnot(!is.na(tad_coexprDT$coexpr))
      stopifnot(is.numeric(tad_coexprDT$coexpr))
      stopifnot(tad_genes %in% tad_coexprDT$gene1 | tad_genes %in% tad_coexprDT$gene2 )
      
      within_coexprDT <- tad_coexprDT[tad_coexprDT$gene1 %in% tad_genes & tad_coexprDT$gene2 %in% tad_genes,]
      between_coexprDT <- tad_coexprDT[ ! (tad_coexprDT$gene1 %in% tad_genes & tad_coexprDT$gene2 %in% tad_genes),]
      
      stopifnot(nrow(within_coexprDT) + nrow(between_coexprDT) == nrow(tad_coexprDT))
      
      list(withinCoexpr = mean(within_coexprDT$coexpr),
           betweenCoexpr = mean(between_coexprDT$coexpr)
           )
    } # end-foreach iterating over TADs
    names(within_between_coexpr_data) <- pipeline_regionList
    within_between_coexpr_data
  } # end-foreach iterating over dataset
  names(allData_within_between_coexpr) <- coexprFilesOutNames
  
  outFile <- file.path(outFolder, "allData_within_between_coexpr.Rdata")
  save(allData_within_between_coexpr, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))
  
} else { # end if-buildData
  outFile <- file.path(outFolder, "allData_within_between_coexpr.Rdata")
  allData_within_between_coexpr <- eval(parse(text = load(outFile)))
}



# ######################################################################################
# ######################################################################################
# ######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))





