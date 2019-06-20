
# Rscript meanCorr_empPval.R

require(foreach)
require(doMC)
registerDoMC(40)

source("utils_fct.R")


script_name <- "meanCorr_empPval.R"
cat("... start ", script_name, "\n")
startTime <- Sys.time()


outFolder <- paste0("MEANCORR_EMPPVAL")
dir.create(outFolder, recursive = TRUE)

acrossFile <- "CORR_AROUND_BD/all_ds_corrData.Rdata"
stopifnot(file.exists(acrossFile))

script4_name <- "4_runMeanTADCorr"

all_ds_corr_across <- eval(parse(text = load(acrossFile)))

all_ds <- names(all_ds_corr_across)

minTADnGenes <- 3

all_ngenes <- 1:10

nameVec <- paste0("nGenes",all_ngenes)

# 
# dist_across_meanCorr <- foreach(ngene = nameVec) %do%{
#   tmp <- unlist(lapply(across_meanCorr, function(x) x[ngene]))
#   names(tmp) <- names(across_meanCorr)
#   tmp
# }
# names(dist_across_meanCorr) <- nameVec

dist_across_meanCorr <- foreach(ngene = nameVec) %do%{
 lapply(all_ds_corr_across, function(nest_list) lapply(nest_list, function(x) x[ngene]))
}
names(dist_across_meanCorr) <- nameVec

ds=all_ds[1]
ds=all_ds[2]

all_ds_empPvals <- foreach(ds = all_ds) %do% {
  
  cat("... start dataset: ", ds, "\n")

  hicds <- dirname(ds)
  stopifnot(dir.exists(hicds))
  
  ds_outFolder <- file.path("PIPELINE", "OUTPUT_FOLDER", ds)
  stopifnot(dir.exists(ds_outFolder))
  obs_meanCorrFile <- file.path(ds_outFolder, script4_name, "all_meanCorr_TAD.Rdata")
  stopifnot(file.exists(obs_meanCorrFile))
  obs_meanCorr <- eval(parse(text=load(obs_meanCorrFile)))
  head(obs_meanCorr)
  
  ### KEEP ONLY THE TADs USED IN THE PIPELINE
  script0_name <- "0_prepGeneData"
  stopifnot(dir.exists(file.path(ds_outFolder, script0_name)))
  tadListFile <- file.path(ds_outFolder, script0_name, "pipeline_regionList.Rdata")
  stopifnot(file.exists(tadListFile))
  pipeline_tadList <- eval(parse(text = load(tadListFile))) # not adjusted

  ### RETRIEVE THE GENES USED IN THE PIPELINE - script0
  stopifnot(dir.exists(file.path(ds_outFolder, script0_name)))
  geneListFile <- file.path(ds_outFolder, script0_name, "pipeline_geneList.Rdata")
  stopifnot(file.exists(geneListFile))
  pipeline_geneList <- eval(parse(text = load(geneListFile))) # not adjusted

  g2tFile <- file.path(hicds, "genes2tad", "all_genes_positions.txt")
  stopifnot(file.exists(g2tFile))
  g2t_DT <- read.delim(g2tFile, header=F, col.names = c("entrezID",  "chromo", "start", "end", "region"), stringsAsFactors = FALSE)
  g2t_DT$entrezID <- as.character(g2t_DT$entrezID)
  
  stopifnot(pipeline_geneList %in% g2t_DT$entrezID)
  
  all_regions <- names(obs_meanCorr)
  reg = all_regions[1]
  reg = all_regions[125]
  
  all_reg_empPvals <- foreach(reg = all_regions) %dopar% {
    
    cat("...... start dataset: ", ds, " - region ", reg, "\n")
    
    
    nGenes <- sum(g2t_DT$entrezID %in% pipeline_geneList & g2t_DT$region == reg)
    stopifnot(nGenes > 0)
    stopifnot(nGenes >= minTADnGenes)

    if(nGenes > 20) {
      # compute the empirical pvalue
      return(list(empPval_allDS = "too_much_genes",
           empPvall_currDS = "too_much_genes"
      ))
    }
        
    if(nGenes %% 2 == 0) {
      stopifnot( paste0("nGenes", nGenes/2) %in% names(dist_across_meanCorr))
      stopifnot(ds %in% names(dist_across_meanCorr[[paste0("nGenes", nGenes/2)]]) )
      
      all_values_allDS <- unlist(dist_across_meanCorr[[paste0("nGenes", nGenes/2)]])
      all_values_currDS <- unlist(dist_across_meanCorr[[paste0("nGenes", nGenes/2)]][[paste0(ds)]])
      
    } else {
      stopifnot( paste0("nGenes", floor(nGenes/2)) %in% names(dist_across_meanCorr))
      stopifnot(ds %in% names(dist_across_meanCorr[[paste0("nGenes", floor(nGenes/2))]]) )
      
      stopifnot( paste0("nGenes", ceiling(nGenes/2)) %in% names(dist_across_meanCorr))
      stopifnot(ds %in% names(dist_across_meanCorr[[paste0("nGenes", ceiling(nGenes/2))]]) )
      
      lower_values_allDS <- unlist(dist_across_meanCorr[[paste0("nGenes", floor(nGenes/2))]])
      upper_values_allDS <- unlist(dist_across_meanCorr[[paste0("nGenes", ceiling(nGenes/2))]])
      all_values_allDS <- c(lower_values_allDS, upper_values_allDS)

      lower_values_currDS <- unlist(dist_across_meanCorr[[paste0("nGenes", floor(nGenes/2))]][[paste0(ds)]])
      upper_values_currDS <- unlist(dist_across_meanCorr[[paste0("nGenes", ceiling(nGenes/2))]][[paste0(ds)]])
      all_values_currDS <- c(lower_values_currDS, upper_values_currDS)
    }
    stopifnot(length(all_values_currDS) > 0)  
    stopifnot(length(all_values_allDS) > 0)  
    
    all_values_allDS <- na.omit(all_values_allDS)
    stopifnot(is.numeric(all_values_allDS))
    
    all_values_currDS <- na.omit(all_values_currDS)
    stopifnot(is.numeric(all_values_currDS))
    
    obs_corrVal <- as.numeric(obs_meanCorr[reg])
    stopifnot(!is.na(obs_corrVal))
    
    # compute the empirical pvalue
    list(empPval_allDS = (sum(all_values_allDS >= obs_corrVal, na.rm=TRUE )+1)/(length(all_values_allDS) + 1),
         empPvall_currDS = (sum(all_values_currDS >= obs_corrVal, na.rm=TRUE )+1)/(length(all_values_currDS) + 1)
    )
  } # end-foreach-iterating over regions
  names(all_reg_empPvals) <- all_regions
  all_reg_empPvals
  
} # end-for-iterating over dataset
names(all_ds_empPvals) <- all_ds

outFile <- file.path(outFolder, "all_ds_empPvals.Rdata")
save(all_ds_empPvals, file = outFile)
cat(paste0("... written: ", outFile, "\n"))




txt <- paste0(startTime, "\n", Sys.time(), "\n")
cat(paste0(txt))
cat(paste0("*** DONE: ", script_name, "\n"))












