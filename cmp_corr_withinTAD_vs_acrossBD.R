source("utils_fct.R")


script_name <- "cmp_corr_withinTAD_vs_acrossBD.R"
cat("... start ", script_name, "\n")
startTime <- Sys.time()


acrossFile <- "CORR_AROUND_BD/all_ds_corrData.Rdata"
stopifnot(file.exists(acrossFile))

script4_name <- "4_runMeanTADCorr"

all_ds_corr_across <- eval(parse(text = load(acrossFile)))

all_ds <- names(all_ds_corr_across)

minTADnGenes <- 3

all_ngenes <- 1:10

nameVec <- paste0("nGenes",all_ngenes)

ds=all_ds[1]
for(ds in all_ds) {
  
  cat("... start dataset: ", ds, "\n")

  hicds <- dirname(ds)
  stopifnot(dir.exists(hicds))
  
  across_meanCorr <- all_ds_corr_across[[ds]]
  nameVec <- names(across_meanCorr[[1]])
  
  dist_across_meanCorr <- foreach(ngene = nameVec) %do%{
    
    tmp <- unlist(lapply(across_meanCorr, function(x) x[ngene]))
    names(tmp) <- names(across_meanCorr)
    tmp
    
  }
  names(dist_across_meanCorr) <- nameVec
  # prep distribution
  plot_multiDens_argList(dist_across_meanCorr)
  plot_multiDens_argList(dist_across_meanCorr[1:5])
  plot_multiDens_argList(dist_across_meanCorr[5:10])
  
    
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
  
  all_reg_empPvals <- foreach(reg = all_regions) %dopar% {
    nGenes <- sum(g2t_DT$entrezID %in% pipeline_geneList & g2t_DT$region == reg)
    stopifnot(nGenes > 0)
    stopifnot(nGenes >= minTADnGenes)
    
    
    if(nGenes %% 2 == 0) {
      all_values <- dist_across_meanCorr[[paste0("nGenes", nGenes/2)]]
      
      
      tmp_list <- list(na.omit(dist_across_meanCorr[[paste0("nGenes", nGenes/2)]]))
      names(tmp_list) <- c(paste0("nGenes", (nGenes/2)))
      plot_multiDens_argList(tmp_list)
      abline(v = obs_corrVal, col="red")
      
      
    } else {
      lower_values <- dist_across_meanCorr[[paste0("nGenes", floor(nGenes/2))]]
      upper_values <- dist_across_meanCorr[[paste0("nGenes", ceiling(nGenes/2))]]
      all_values <- c(lower_values, upper_values)
      
      tmp_list <- list(dist_across_meanCorr[[paste0("nGenes", floor(nGenes/2))]],
                        dist_across_meanCorr[[paste0("nGenes", ceiling(nGenes/2))]])
      names(tmp_list) <- c(paste0("nGenes", floor(nGenes/2)), paste0("nGenes", ceiling(nGenes/2)))
      
      plot_multiDens_argList(tmp_list)
      abline(v = obs_corrVal, col="red")
      
      
    }
    stopifnot(length(all_values) > 0)    
    
    all_values <- na.omit(all_values)
    
    obs_corrVal <- as.numeric(obs_meanCorr[reg])
    stopifnot(!is.na(obs_corrVal))
    
    # compute the empirical pvalue
    (sum(all_values >= obs_corrVal, na.rm=TRUE )+1)/(length(all_values) + 1)
  } # end-foreach-iterating over regions
  
  
} # end-for-iterating over dataset



tadSizes <- foreach(reg = pipeline_regionList, .combine='c') %dopar% {
  sum(g2t_DT$entrezID %in% pipeline_geneList & g2t_DT$region == reg)
}
names(tadSizes) <- pipeline_regionList
stopifnot(tadSizes >= minTADnGenes)

nGenes = 6

for(nGenes in minTADnGenes:(2*max(all_ngenes))) {
  
  
  obs_tads <- names(tadSizes)[tadSizes == nGenes]
  obs_corrs <- as.numeric(obs_meanCorr[names(obs_meanCorr) %in% obs_tads])
  stopifnot(!is.na(obs_corrs))
  

    
    if(nGenes %% 2 == 0) {
      all_values <- dist_across_meanCorr[[paste0("nGenes", nGenes/2)]]
      
      tmp_list <- list(na.omit(dist_across_meanCorr[[paste0("nGenes", nGenes/2)]]),
                       obs_corrs)
      names(tmp_list) <- c(paste0("nGenes", (nGenes/2)), "obs. values")
      plot_multiDens_argList(tmp_list)
      
      
    } else {
      lower_values <- dist_across_meanCorr[[paste0("nGenes", floor(nGenes/2))]]
      upper_values <- dist_across_meanCorr[[paste0("nGenes", ceiling(nGenes/2))]]
      all_values <- c(lower_values, upper_values)
      
      tmp_list <- list(dist_across_meanCorr[[paste0("nGenes", floor(nGenes/2))]],
                       dist_across_meanCorr[[paste0("nGenes", ceiling(nGenes/2))]],
                       obs_corrs)
      names(tmp_list) <- c(paste0("nGenes", floor(nGenes/2)), 
                           paste0("nGenes", ceiling(nGenes/2)),
                           "obs. values")
      
      plot_multiDens_argList(tmp_list)

            
      
    }
  stopifnot(length(all_values) > 0)    
  
  all_values <- na.omit(all_values)
  
  
  
  
}




txt <- paste0(startTime, "\n", Sys.time(), "\n")
cat(paste0(txt))
cat(paste0("*** DONE: ", script_name, "\n"))












