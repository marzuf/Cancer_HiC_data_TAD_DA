# Rscript corr_within_TAD.R

# for a given number of genes
# mean corr within TAD

script_name <- "corr_within_TAD.R"

cat("... start ", script_name, "\n")

startTime <- Sys.time()

SSHFS <- FALSE
nCpu <- ifelse(SSHFS, 2, 40)

nCpu <- 1

suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
registerDoMC(nCpu)

mywd <- ifelse(SSHFS, "/media/electron/mnt/etemp/marie/Cancer_HiC_data_TAD_DA", "")

outFolder <- file.path("CORR_WITHIN_TAD")
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

# all_gene_nbrs <- 1:10
all_gene_nbrs <- 1:3

all_ds_corrData <- list()

# all_hicexpr_ds=all_hicexpr_ds[1]
ds=all_hicexpr_ds[1]
all_hicexpr_ds=all_hicexpr_ds[1:2]
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
  
  
  ### RETRIEVE THE GENES USED IN THE PIPELINE - script0
  script0_name <- "0_prepGeneData"
  stopifnot(dir.exists(file.path(dsPipOutDir, script0_name)))
  geneListFile <- file.path(dsPipOutDir, script0_name, "pipeline_geneList.Rdata")
  stopifnot(file.exists(geneListFile))
  pipeline_geneList <- eval(parse(text = load(geneListFile))) # not adjusted
  stopifnot(pipeline_geneList %in% g2t_DT$entrezID)
  
  g2t_DT <- g2t_DT[g2t_DT$entrezID %in% as.character(pipeline_geneList),]
  
  norm_rnaseqDT <- eval(parse(text = load(file.path(dsPipOutDir, script0_name, "rna_qqnorm_rnaseqDT.Rdata"))))
  stopifnot(names(pipeline_geneList) %in% rownames(norm_rnaseqDT))
  # stopifnot(pipeline_geneList %in% rownames(norm_rnaseqDT)) -> FALSE
  norm_rnaseqDT <- norm_rnaseqDT[names(pipeline_geneList),]    
  stopifnot(rownames(norm_rnaseqDT) == names(pipeline_geneList))
  

  all_regions <- unique(tadpos_DT$region)
  
  stopifnot(sort(all_regions) == sort(pipeline_regionList))

  
  all_regions=all_regions[1:2]
  
  nTADs <- length(all_regions)
  reg = all_regions[1]
  all_meanCorr_TAD <- foreach(reg=all_regions) %dopar% {

    reg_genes <- g2t_DT$entrezID[g2t_DT$region == reg]
    stopifnot(length(reg_genes) > 0)
    
    stopifnot(reg_genes %in% pipeline_geneList)
    
    
    tad_all_nGenes_corr <- rep(NA, length(all_gene_nbrs))
    names(tad_all_nGenes_corr) <- as.character(2*all_gene_nbrs)
    
    
    nGenes = 1
    for(nGenes in all_gene_nbrs) {
      
      n_genes <- nGenes*2 # because I chose nGenes right + nGenes left
      
      cat("......... start # genes = ", n_genes, "\n")
      
      
      if( n_genes > length(reg_genes)) {
        tad_all_nGenes_corr[as.character(n_genes)] <- NA
        next
      }
      
      
      
      all_cmbs_DT <- combn(reg_genes, n_genes)
      
      stopifnot(nrow(all_cmbs_DT) == n_genes )
      
      nCmbs <- ncol(all_cmbs_DT)
      
      n <- length(reg_genes)
      k <- n_genes
      
      stopifnot( factorial(n)/ (factorial(k)*factorial(n-k)) == nCmbs)
      
      
       
      i_cmb=1
      all_cmbs <- foreach(i_cmb = seq_len(nCmbs), .combine='c') %do% {
        
        cat("............ start cmb idx = ", i_cmb, "\n")
        
        
        cmb_reg_genes <- all_cmbs_DT[,i_cmb]
        stopifnot(length(cmb_reg_genes) == n_genes)
        stopifnot(cmb_reg_genes %in% reg_genes)
        stopifnot(cmb_reg_genes %in% pipeline_geneList)
        stopifnot(cmb_reg_genes %in% g2t_DT$entrezID)
        rnaToKeep <- names(pipeline_geneList)[pipeline_geneList %in% cmb_reg_genes]
        stopifnot(rnaToKeep %in% rownames(rna_qqnorm_rnaseqDT))
        sub_expr_DT <- rna_qqnorm_rnaseqDT[rnaToKeep,,drop=FALSE]
        corMatrixWithinTAD <- cor(t(sub_expr_DT), method=corMet)
        stopifnot(dim(corMatrixWithinTAD) == nGenes*2)
        
        meanCorMatrixWithinTAD <- mean(corMatrixWithinTAD[lower.tri(corMatrixWithinTAD, diag = withDiago)], na.rm=TRUE)
        
        meanCorMatrixWithinTAD
        
        
        
      } # end-for each combination of nGenes genes possible within current TAD 
      stopifnot(length(all_cmbs) == ncol(all_cmbs_DT))
      stopifnot(!is.na(all_cmbs))
      nGenes_meanCorr <- mean(all_cmbs, na.rm = TRUE)
      tad_all_nGenes_corr[as.character(n_genes)] <- nGenes_meanCorr
      
    } # end-for iterating over number of genes
    stopifnot(length(tad_all_nGenes_corr) == length(all_gene_nbrs))
    #stopifnot(!is.na(bd_all_nGenes_corr)) # not true if not enough genes on the right or on the right
    
    names(tad_all_nGenes_corr) <- paste0("nGenes", names(tad_all_nGenes_corr))
    tad_all_nGenes_corr

          
  } # end-foreach-iterating over the TADs
  
  stopifnot(length(all_meanCorr_TAD) == nTADs)
  names(all_meanCorr_TAD) <- paste0("TAD", 1:nTADs)
  all_ds_corrData[[paste0(ds)]] <- all_meanCorr_TAD
  
  
  
  
} # end-for iterating over datasets

outFile <- file.path(outFolder, "all_ds_corrData.Rdata")
save(all_ds_corrData, file=outFile)

cat(paste0("... written: ", outFile, "\n"))

txt <- paste0(startTime, "\n", Sys.time(), "\n")
cat(paste0(txt))
cat(paste0("*** DONE: ", script_name, "\n"))







