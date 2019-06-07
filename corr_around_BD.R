# Rscript corr_around_BD.R

# for a given number of genes
# go at each BD location
# select "n" genes left and right of the BD
# to get an empirical distribution of the correlation
# between expression of genes separated by the boundary

script_name <- "corr_around_BD.R"

startTime <- Sys.time()

SSHFS <- FALSE
nCpu <- ifelse(SSHFS, 2, 40)

suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
registerDoMC(nCpu)

mywd <- ifelse(SSHFS, "/media/electron/mnt/etemp/marie/Cancer_HiC_data_TAD_DA", "")

outFolder <- file.path("CORR_AROUND_BD")
dir.create(outFolder)

pipOutFolder <- file.path("PIPELINE", "OUTPUT_FOLDER")

stopifnot(dir.exists(pipOutFolder))

all_hicexpr_ds <- unname(unlist(sapply(list.files(pipOutFolder, full.names = TRUE), function(x) file.path(basename(x),  list.files(x)))))
stopifnot(length(all_hicexpr_ds) > 0)
stopifnot(dir.exists(file.path(pipOutFolder, all_hicexpr_ds)))

ds="ENCSR862OGI_RPMI-7951_40kb/TCGAskcm_lowInf_highInf"
ds=all_hicexpr_ds[1]

corMet <- "pearson"
withDiago <- FALSE

all_gene_nbrs <- 1:10

all_ds_corrData <- list()

# all_hicexpr_ds=all_hicexpr_ds[1]
for(ds in all_hicexpr_ds) {
  hicds <- file.path(dirname(ds))
  exprds <- basename(ds)
  stopifnot(dir.exists(hicds))
  
  cat("... start dataset: ", exprds, "\n")
  
  dsPipOutDir <- file.path(pipOutFolder, ds)
  stopifnot(dir.exists(dsPipOutDir))
  
  ### RETRIEVE THE GENE2TAD ASSIGNMENT
  g2tFile <- file.path(hicds, "genes2tad", "all_genes_positions.txt")
  stopifnot(file.exists(g2tFile))
  g2t_DT <- read.delim(g2tFile, header=F, col.names = c("entrezID",  "chromo", "start", "end", "region"), stringsAsFactors = FALSE)
  g2t_DT$entrezID <- as.character(g2t_DT$entrezID)
  
  ### RETRIEVE THE TAD POSITIONS
  tadposFile <- file.path(hicds, "genes2tad", "all_assigned_regions.txt")
  stopifnot(file.exists(tadposFile))
  tadpos_DT <- read.delim(tadposFile, header=F, col.names=c("chromo", "region", "start", "end"), stringsAsFactors = FALSE)
  stopifnot(is.numeric(tadpos_DT$start))
  stopifnot(is.numeric(tadpos_DT$end))
  tadpos_DT <- tadpos_DT[grepl("_TAD", tadpos_DT$region),,drop=FALSE] 
  
  
  ### KEEP ONLY THE TADs USED IN THE PIPELINE
  script0_name <- "0_prepGeneData"
  stopifnot(dir.exists(file.path(dsPipOutDir, script0_name)))
  tadListFile <- file.path(dsPipOutDir, script0_name, "pipeline_regionList.Rdata")
  stopifnot(file.exists(tadListFile))
  pipeline_tadList <- eval(parse(text = load(tadListFile))) # not adjusted
  stopifnot(pipeline_tadList %in% tadpos_DT$region)
  tadpos_DT <- tadpos_DT[tadpos_DT$region %in% pipeline_tadList,]
  
  
  bdpos_DT1 <- tadpos_DT[, c("chromo", "start")]
  bdpos_DT1$start <- bdpos_DT1$start - 1
  colnames(bdpos_DT1) <- c("chromo", "BDpos")
  bdpos_DT2 <- tadpos_DT[, c("chromo", "end")]
  colnames(bdpos_DT2) <- c("chromo", "BDpos")
  bdpos_DT <- rbind(bdpos_DT1, bdpos_DT2)
  stopifnot(nrow(bdpos_DT) == 2*nrow(tadpos_DT))
  bdpos_DT <- unique(bdpos_DT)
  nPos <- nrow(bdpos_DT)
  
  ### RETRIEVE THE GENES USED IN THE PIPELINE - script0
  script0_name <- "0_prepGeneData"
  stopifnot(dir.exists(file.path(dsPipOutDir, script0_name)))
  geneListFile <- file.path(dsPipOutDir, script0_name, "pipeline_geneList.Rdata")
  stopifnot(file.exists(geneListFile))
  pipeline_geneList <- eval(parse(text = load(geneListFile))) # not adjusted
  stopifnot(pipeline_geneList %in% g2t_DT$entrezID)
  
  stopifnot(pipeline_geneList %in% g2t_DT$entrezID)
  # stopifnot(names(pipeline_geneList) %in% g2t_DT$entrezID) -> FALSE
  g2t_DT <- g2t_DT[as.character(g2t_DT$entrezID) %in% as.character(pipeline_geneList),,drop=FALSE]
  stopifnot(length(pipeline_geneList) == nrow(g2t_DT))
  stopifnot(g2t_DT$entrezID %in% pipeline_geneList)
  g2t_DT$chromo <- as.character(g2t_DT$chromo)
  
  stopifnot(g2t_DT$entrezID %in% pipeline_geneList)
  stopifnot(grepl("TAD", g2t_DT$region))
  
  norm_rnaseqDT <- eval(parse(text = load(file.path(dsPipOutDir, script0_name, "rna_qqnorm_rnaseqDT.Rdata"))))
  stopifnot(names(pipeline_geneList) %in% rownames(norm_rnaseqDT))
  # stopifnot(pipeline_geneList %in% rownames(norm_rnaseqDT)) -> FALSE
  norm_rnaseqDT <- norm_rnaseqDT[names(pipeline_geneList),]    
  stopifnot(rownames(norm_rnaseqDT) == names(pipeline_geneList))
  
  # nPos=10
  allBD_nGenes_corr <- foreach(i_bd = seq_len(nPos)) %dopar% {
  # allBD_nGenes_corr <- foreach(i_bd = seq_len(nPos)) %do% {
    
    cat("...... start BD pos. : \t", i_bd, "/", nPos, "\n")
    
    
    curr_chromo <- as.character(bdpos_DT$chromo[i_bd])
  
    curr_pos <- bdpos_DT$BDpos[i_bd]
    stopifnot(is.numeric(curr_pos))
    stopifnot(length(curr_pos) == 1)
    
    # !!! EXTRACT GENES BASED ON START POSITION RELATIVE TO BD
    # !!! SMALLER THAN / GREATER *OR EQUAL* THAN BD POSITION (smaller not equal otherwise genes could come twice)
    
    curr_g2t <- g2t_DT[g2t_DT$chromo == curr_chromo,,drop=FALSE]
    stopifnot(nrow(curr_g2t) > 0)
    
    stopifnot(is.numeric(curr_g2t$start), is.numeric(curr_g2t$end))
    curr_g2t <- curr_g2t[order(curr_g2t$start, curr_g2t$end),,drop=FALSE]
    
    
    curr_g2t$posDiff <- curr_pos - curr_g2t$start
    
    curr_genesLeftDT <- curr_g2t[curr_g2t$start < curr_pos,,drop=FALSE] 
    stopifnot(nrow(curr_genesLeftDT) == sum(curr_g2t$posDiff > 0))
    
    curr_genesRightDT <- curr_g2t[curr_g2t$start >= curr_pos,,drop=FALSE] 
    stopifnot(nrow(curr_genesRightDT) == sum(curr_g2t$posDiff <= 0))
    
    nrow(curr_genesRightDT)
    nrow(curr_genesLeftDT)
    
    if(nrow(curr_genesRightDT) == 0 | nrow(curr_genesLeftDT) == 0) {
      return(NA)
    }

    stopifnot(curr_genesRightDT$entrezID %in% pipeline_geneList)   
    stopifnot(curr_genesLeftDT$entrezID %in% pipeline_geneList)   
    
    # on the right: negative, sort from biggest to smallest (decreasing)
    # on the left: positive, sort from smallest to biggest (increasing)
    curr_genesRightDT <- curr_genesRightDT[order(curr_genesRightDT$posDiff, decreasing = TRUE),,drop=FALSE]
    curr_genesLeftDT <- curr_genesLeftDT[order(curr_genesLeftDT$posDiff, decreasing = FALSE),,drop=FALSE]
    
    stopifnot(is.numeric(curr_genesLeftDT$posDiff))
    stopifnot(is.numeric(curr_genesRightDT$posDiff))
    
    # > all(pipeline_geneList %in% rownames(rna_qqnorm_rnaseqDT))
    # [1] FALSE
    # > all(names(pipeline_geneList) %in% rownames(rna_qqnorm_rnaseqDT))
    # [1] TRUE
    
    # stopifnot(as.character(curr_genesRightDT$entrezID) %in% as.character(rownames(rna_qqnorm_rnaseqDT)))
    # stopifnot(as.character(curr_genesLeftDT$entrezID) %in% as.character(rownames(rna_qqnorm_rnaseqDT)))
    
    stopifnot(names(pipeline_geneList)[pipeline_geneList %in% as.character(curr_genesRightDT$entrezID)] %in%
                as.character(rownames(rna_qqnorm_rnaseqDT)))
    stopifnot(names(pipeline_geneList)[pipeline_geneList %in% as.character(curr_genesLeftDT$entrezID)] %in% 
                as.character(rownames(rna_qqnorm_rnaseqDT)))
    
    #all_gene_nbrs=1:10
    
    bd_all_nGenes_corr <- rep(NA, length(all_gene_nbrs))
    names(bd_all_nGenes_corr) <- as.character(all_gene_nbrs)

    for(nGenes in all_gene_nbrs) {
      
      if( nrow(curr_genesLeftDT) < nGenes | nrow(curr_genesRightDT) < nGenes ) {
        bd_all_nGenes_corr[as.character(nGenes)] <- NA
        next
      }
      
      cat("......... start # genes = ", nGenes, "\n")
      
      # select nGenes numbers from each side of each boundary
      # > all(curr_g2t$entrezID %in% names(pipeline_geneList))
      # [1] FALSE
      # > all(curr_g2t$entrezID %in% (pipeline_geneList))
      # [1] TRUE
      # 
      # entrezID <-> geneList
      # names(geneList) <-> rownames rnaseq
      
      stopifnot(curr_genesLeftDT$entrezID %in% pipeline_geneList)
      stopifnot(curr_genesRightDT$entrezID %in% pipeline_geneList)
      
      nGenesLeft <- curr_genesLeftDT$entrezID[1:nGenes]
      stopifnot(nGenesLeft %in% pipeline_geneList)
      # stopifnot(nGenesLeft %in% names(pipeline_geneList))
      nGenesLeft_rnaID <- names(pipeline_geneList)[pipeline_geneList %in% nGenesLeft]
      stopifnot(nGenesLeft_rnaID %in% rownames(rna_qqnorm_rnaseqDT))
      
      stopifnot(!is.na(nGenesLeft))
      stopifnot(is.character(nGenesLeft))
      
      nGenesRight <- curr_genesRightDT$entrezID[1:nGenes]
      stopifnot(nGenesRight %in% pipeline_geneList)
      # stopifnot(nGenesRight %in% names(pipeline_geneList))
      nGenesRight_rnaID <- names(pipeline_geneList)[pipeline_geneList %in% nGenesRight]
      stopifnot(nGenesRight_rnaID %in% rownames(rna_qqnorm_rnaseqDT))
      
      stopifnot(!nGenesRight %in% nGenesLeft)
      
      stopifnot(!is.na(nGenesRight))
      stopifnot(is.character(nGenesRight))
      
      stopifnot(length(nGenesLeft) == nGenes)
      stopifnot(length(nGenesRight) == nGenes)
    
      ### OR PUT AS A SINGLE VECTOR ???
      exprLeft <- rna_qqnorm_rnaseqDT[as.character(nGenesLeft_rnaID),,drop=FALSE]
      exprRight <- rna_qqnorm_rnaseqDT[as.character(nGenesRight_rnaID),,drop=FALSE]
      stopifnot(nrow(exprLeft) == nrow(exprRight))
      stopifnot(dim(exprLeft) == dim(exprRight))
      # cor(as.numeric(exprLeft), as.numeric(exprRight))
          
      exprRightLeft <- rna_qqnorm_rnaseqDT[c(as.character(nGenesLeft_rnaID), as.character(nGenesRight_rnaID)),,drop=FALSE]
      stopifnot(nrow(exprRightLeft) == nrow(exprLeft) + nrow(exprRight))
    
      corMatrixRightLeft <- cor(t(exprRightLeft), method=corMet)
      stopifnot(dim(corMatrixRightLeft) == nGenes*2)
      
      meanCorr_rightLeft <- mean(corMatrixRightLeft[lower.tri(corMatrixRightLeft, diag = withDiago)], na.rm=TRUE)
      
      bd_all_nGenes_corr[as.character(nGenes)] <- meanCorr_rightLeft # should match the name used when create vector full of NA
      
    } # end for-iterating over gene #
    stopifnot(length(bd_all_nGenes_corr) == length(all_gene_nbrs))
    #stopifnot(!is.na(bd_all_nGenes_corr)) # not true if not enough genes on the right or on the right
    
    names(bd_all_nGenes_corr) <- paste0("nGenes", names(bd_all_nGenes_corr))
    bd_all_nGenes_corr
    
  } # end foreach-iterating over boundary positions
  stopifnot(length(allBD_nGenes_corr) == nPos)
  names(allBD_nGenes_corr) <- paste0("boundary", 1:nPos)
  all_ds_corrData[[paste0(ds)]] <- allBD_nGenes_corr
} # end for iterating over ds

  
outFile <- file.path(outFolder, "all_ds_corrData.Rdata")
save(all_ds_corrData, file=outFile)

cat(paste0("... written: ", outFile, "\n"))

txt <- paste0(startTime, "\n", Sys.time(), "\n")
cat(paste0(txt))
cat(paste0("*** DONE: ", script_name, "\n"))











