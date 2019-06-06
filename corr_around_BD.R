# for a given number of genes
# go at each BD location
# select "n" genes left and right of the BD
# to get an empirical distribution of the correlation
# between expression of genes separated by the boundary

mywd <- ""
mywd <- "/media/electron/mnt/etemp/marie/Cancer_HiC_data_TAD_DA"

pipOutFolder <- file.path(mywd, "PIPELINE", "OUTPUT_FOLDER")

all_hicexpr_ds <- unname(unlist(sapply(list.files(pipOutFolder, full.names = TRUE), function(x) file.path(basename(x),  list.files(x)))))
stopifnot(dir.exists(file.path(pipOutFolder, all_hicexpr_ds)))

ds="ENCSR862OGI_RPMI-7951_40kb/TCGAskcm_lowInf_highInf"

all_gene_nbrs <- 1:10

for(ds in all_hicexpr_ds) {
  hicds <- file.path(mywd, dirname(ds))
  exprds <- basename(ds)
  stopifnot(dir.exists(hicds))
  
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
  tadpos_DT <- tadpos_DT[grepl("_TAD", tadpos_DT$region),] 
  
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
  g2t_DT <- g2t_DT[g2t_DT$entrezID %in% pipeline_geneList,]
  stopifnot(length(pipeline_geneList) == nrow(g2t_DT))
  g2t_DT$chromo <- as.character(g2t_DT$chromo)
  
  stopifnot(grepl("TAD", g2t_DT$region))
  
  for(i_bd in seq_len(nPos)) {
    
    cat("... start BD pos. : \t", i_bd, "/", nPos, "\n")
    
    
    curr_chromo <- as.character(bdpos_DT$chromo[i_bd])
  
    curr_pos <- bdpos_DT$BDpos[i_bd]
    stopifnot(is.numeric(curr_pos))
    stopifnot(length(curr_pos) == 1)
    
    # !!! EXTRACT GENES BASED ON START POSITION RELATIVE TO BD
    # !!! SMALLER THAN / GREATER *OR EQUAL* THAN BD POSITION (smaller not equal otherwise genes could come twice)
    
    curr_g2t <- g2t_DT[g2t_DT$chromo == curr_chromo,]
    stopifnot(nrow(curr_g2t) > 0)
    
    stopifnot(is.numeric(curr_g2t$start), is.numeric(curr_g2t$end))
    curr_g2t <- curr_g2t[order(curr_g2t$start, curr_g2t$end),]
    
    
    curr_g2t$posDiff <- curr_pos - curr_g2t$start
    
    curr_genesLeftDT <- curr_g2t[curr_g2t$start < curr_pos,] 
    stopifnot(nrow(curr_genesLeftDT) == sum(curr_g2t$posDiff > 0))
    
    
    curr_genesRightDT <- curr_g2t[curr_g2t$start >= curr_pos,] 
    stopifnot(nrow(curr_genesRightDT) == sum(curr_g2t$posDiff <= 0))
    
    
    if(nrow(curr_genesRightDT) == 0 | nrow(curr_genesLeftDT) == 0) {
      next
    }
    
    
    curr_genesRightDT <- curr_genesRightDT[order(curr_genesRightDT$posDiff),]
    curr_genesLeftDT <- curr_genesLeftDT[order(curr_genesLeftDT$posDiff),]
    
    stopifnot(is.numeric(curr_genesLeftDT$posDiff))
    stopifnot(is.numeric(curr_genesRightDT$posDiff))
    
    
    # take pairwise correlation, or concatenate the vectors ???
    
    nGenes <- 2
    for(nGenes in all_gene_nbrs) {
      
      # select nGenes numbers from each side of each boundary
      
      
      nGenesLeft <- curr_genesLeftDT$entrezID[1:nGenes]
      
      stopifnot(!is.na(nGenesLeft))
      stopifnot(is.character(nGenesLeft))
      
      nGenesRight <- curr_genesRightDT$entrezID[1:nGenes]
      
      stopifnot(!is.na(nGenesRight))
      stopifnot(is.character(nGenesRight))
      
      
      
      
      
    }
    
    
  }
  
  
  ### ITERATE OVER THE CHROMOSOMES !!!
  all_chr <- unique(as.character(tadpos_DT$chromo))
  
  chromo <-  "chr1"
  for(chromo in all_chr) {
    
    chromoTADpos_DT <- tadpos_DT[tadpos_DT$chromo == chromo,]
    stopifnot(nrow(chromoTADpos_DT) > 0)
    
    # remove 1 -> when I take the boundary vector, I don't want "duplicates" (e.g. an end at 50'000 and a start at 50'001)
    chromoTADpos_DT$start <- chromoTADpos_DT$start - 1 
    all_boundaries <- unique(c(chromoTADpos_DT$start, chromoTADpos_DT$end))
    stopifnot(is.numeric(all_boundaries))
    
    # iterate i) over the boundaries; ii) over the nGenes
    
    for(i in seq_along(all_boundaries)) {
      
      
        
    
  }
  
  
  
  
      
    
    
    
    
  }
  
}    
