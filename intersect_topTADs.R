
# Rscript intersect_topTADs.R TCGAkich_norm_kich 3 ENCSR079VIJ_G401_40kb ENCSR401TBQ_Caki2_40kb ENCSR079VIJ_G401ENCSR401TBQ_Caki2_40kb
# Rscript intersect_topTADs.R TCGAluad_mutKRAS_mutEGFR 3 ENCSR444WCZ_A549_40kb NCI-H460_40kb ENCSR444WCZ_A549NCI-H460_40kb
# Rscript intersect_topTADs.R TCGAskcm_wt_mutCTNNB1 3 ENCSR312KHQ_SK-MEL-5_40kb ENCSR862OGI_RPMI-7951_40kb ENCSR312KHQ_SK-MEL-5ENCSR862OGI_RPMI-7951_40kb

# 
# 


startTime <- Sys.time()

require(foreach)
require(IRanges)
require(GenomicRanges)
source("utils_fct.R")

setDir=""
setDir <- "~/media/electron"
setDir="" 

outFolder <- file.path("INTERSECT_topTADs")
dir.create(outFolder, recursive = TRUE)

args <- commandArgs(trailingOnly = TRUE)

all_hicds <- c("ENCSR312KHQ_SK-MEL-5_40kb", "ENCSR862OGI_RPMI-7951_40kb", "ENCSR312KHQ_SK-MEL-5ENCSR862OGI_RPMI-7951_40kb")
exprds <- "TCGAskcm_wt_mutCTNNB1"
topThresh <- 1
topThresh <- 5
all_hicds <- c("ENCSR444WCZ_A549_40kb", "NCI-H460_40kb" ,"ENCSR444WCZ_A549_NCI-H460_40kb", "pipelineConsensus")
exprds <- "TCGAluad_mutKRAS_mutEGFR"
topThresh <- 10
all_hicds <- "ENCSR079VIJ_G401_40kb"
exprds <- "TCGAkich_norm_kich"
topThresh <- 3
stopifnot(length(args) > 2)


all_hicds <- args[3:length(args)]
exprds <- args[1]
topThresh <- as.numeric(args[2])

logFile=""

stopifnot(!is.na(topThresh))
stopifnot(topThresh > 0)

entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
entrezDT <- read.delim(entrezDT_file, header=TRUE, stringsAsFactors = FALSE)
entrezDT$entrezID <- as.character(entrezDT$entrezID)
head(entrezDT)

#PIPELINE/OUTPUT_FOLDER/GSE105318_DLD1_40kb/TCGAcoad_msi_mss//emp_pval_combined.Rdata
script0_name <- "0_prepGeneData"
script11_name <- "11_runEmpPvalCombined"

all_topTADsGenes <- foreach(hicds = all_hicds) %dopar% {
  
  if(hicds == "pipelineConsensus") {
    pipOutDir <- file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom", "OUTPUT_FOLDER",  exprds)
    
    # file with assignment from entrez to all regions
    gene2tadDT_file <- paste0(setDir, 
                              "/mnt/ed4/marie/gene_data_final/consensus_TopDom_covThresh_r0.6_t80000_v0_w-1_final/genes2tad/all_genes_positions.txt")
    TADposDT_file <- paste0(setDir, 
                              "/mnt/ed4/marie/gene_data_final/consensus_TopDom_covThresh_r0.6_t80000_v0_w-1_final/genes2tad/all_assigned_regions.txt") 
    
  } else {
    stopifnot(dir.exists(hicds))
    pipOutDir <- file.path("PIPELINE", "OUTPUT_FOLDER", hicds, exprds)
    gene2tadDT_file <- file.path(hicds, "genes2tad/all_genes_positions.txt")
    TADposDT_file <- file.path(hicds, "genes2tad/all_assigned_regions.txt")
  }
  stopifnot(file.exists(gene2tadDT_file))
  g2tDT <- read.delim(gene2tadDT_file, stringsAsFactors = FALSE, header=F, col.names=c("entrezID", "chromo", "start", "end", "region"))
  
  stopifnot(file.exists(TADposDT_file))
  TADposDT <- read.delim(TADposDT_file, stringsAsFactors = FALSE, header=F, col.names=c("chromo","region", "start", "end"))
  stopifnot(is.numeric(TADposDT$start))
  TADposDT <- TADposDT[grepl("_TAD", TADposDT$region),]
  stopifnot(nrow(TADposDT) > 0)
  
  stopifnot(g2tDT$entrezID %in% entrezDT$entrezID)
  g2t2sDT <- merge(g2tDT[, c("entrezID", "region")], entrezDT[,c("entrezID", "symbol")], by="entrezID", all.x = TRUE, all.y = FALSE )
  head(g2t2sDT)
  stopifnot(!is.na(g2t2sDT))
  g2t2sDT$entrezID <- as.character(g2t2sDT$entrezID)
  g2t2sDT$region <- as.character(g2t2sDT$region)
  g2t2sDT$symbol <- as.character(g2t2sDT$symbol)
  
  # subset the genes used in the pipeline
  stopifnot(dir.exists(pipOutDir))
  stopifnot(dir.exists(file.path(pipOutDir, script0_name)))
  # geneList: values are entrez, names can be ENSEMBL, entrez, etc.
  geneList_file <- file.path(pipOutDir, script0_name, "pipeline_geneList.Rdata")
  stopifnot(file.exists(geneList_file))
  geneList <- eval(parse(text = load(geneList_file)))
  stopifnot(geneList %in% g2tDT$entrezID)
  g2tDT <- g2tDT[g2tDT$entrezID %in% geneList,]
  
  stopifnot(dir.exists(file.path(pipOutDir, script11_name)))
  tad_pvalFile <- file.path(pipOutDir, script11_name, "emp_pval_combined.Rdata")
  stopifnot(file.exists(tad_pvalFile))
  tad_pvals <- eval(parse(text = load(tad_pvalFile)))
  
  adj_tad_pvals <- sort(p.adjust(tad_pvals, method="BH"))
  
  if(topThresh >= 1) {
    topThresh <- min(c(topThresh, length(adj_tad_pvals)))
    pvalThresh <- as.numeric(adj_tad_pvals[topThresh])
    stopifnot(!is.na(pvalThresh))
    stopifnot(pvalThresh <= 1)
  } else {
    pvalThresh <- topThresh
  }
  top_pvals <- adj_tad_pvals[adj_tad_pvals <= pvalThresh]
  stopifnot(!is.na(top_pvals))
  
  if(topThresh >= 1) stopifnot(length(unique(top_pvals)) <= topThresh)
  
  top_tads <- names(top_pvals)
  
  tad_ranks <- rank(top_pvals, ties="min")
  stopifnot(top_tads %in% names(tad_ranks))
  tad_ranks <- tad_ranks[top_tads]
  
  topTADsDT <- data.frame(
    top_tads = top_tads,
    tad_rank = tad_ranks,
    top_pvals = as.numeric(top_pvals),
    stringsAsFactors = FALSE
  )
  head(topTADsDT)
  
  stopifnot(top_tads %in% g2t2sDT$region)
  topGenesDT <- g2t2sDT[g2t2sDT$region %in% top_tads,]
  head(topGenesDT)
  
  colnames(topGenesDT)[colnames(topGenesDT) == "region"] <- "top_tads"
  
  topTADs_genesDT <- merge(topTADsDT, topGenesDT[, c("top_tads", "symbol")], 
                           by = "top_tads", all.x=TRUE, all.y=TRUE)
  
  topTADs_genesDT$top_tads <- as.character(topTADs_genesDT$top_tads)
  topTADs_genesDT$symbol <- as.character(topTADs_genesDT$symbol)
  topTADs_genesDT <- topTADs_genesDT[order(topTADs_genesDT$top_pvals, topTADs_genesDT$top_tads, topTADs_genesDT$symbol),]
  head(topTADs_genesDT)
  
  list(topTADs_genesDT=topTADs_genesDT,
       TADposDT=TADposDT,
       gene2tadDT = g2tDT)
}
names(all_topTADsGenes) <- all_hicds
str(all_topTADsGenes)

########################################################################### LOAD SAVED DATA
# outFolder <- "INTERSECT_topTADs"
# outFile <- file.path(outFolder, "all_topTADsGenes.Rdata")
# save(all_topTADsGenes, file = outFile)
      # outFile = "INTERSECT_topTADs/all_topTADsGenes.Rdata"
      # load(outFile)
      # load("INTERSECT_topTADs/all_topTADsGenes.Rdata")
      # str(all_topTADsGenes)

# names(all_topTADsGenes[[1]])
# [1] "topTADs_genesDT" "TADposDT"        "gene2tadDT" 

i=1
j=2

all_maxIntersectOverlap_matchDT <- foreach(i = seq_along(all_topTADsGenes), .combine="rbind") %dopar% {
  
  ref_ds <- names(all_topTADsGenes)[i]
  
  cat("... start with ref_ds\t=\t", ref_ds, "\n")
  
  signifTADs <- unique(as.character(all_topTADsGenes[[i]][["topTADs_genesDT"]][, "top_tads"]))
  
  signifTADs_posDT <- all_topTADsGenes[[i]][["TADposDT"]]
  signifTADs_posDT$region <- as.character(signifTADs_posDT$region)
  stopifnot(signifTADs %in% signifTADs_posDT$region)
  signifTADs_posDT <- signifTADs_posDT[signifTADs_posDT$region %in% signifTADs,]
  
  ref_g2tDT <- all_topTADsGenes[[i]][["gene2tadDT"]]
  
  query_signifTADs_IR <- IRanges(start = signifTADs_posDT$start, 
                           width = (signifTADs_posDT$end - signifTADs_posDT$start + 1), 
                           names=signifTADs_posDT$region)
  
  
  query_signifTADs_GR <- GRanges(ranges = query_signifTADs_IR,
                                 seqnames=gsub("(chr.+)_TAD.+", "\\1", signifTADs_posDT$region))
  
  
  ref_maxIntersectOverlap_matchDT <- foreach(j = seq_along(all_topTADsGenes)[-i], .combine="rbind") %dopar% {
    
    obj_ds <- names(all_topTADsGenes)[j]
    
    cat("... start with obj_ds\t=\t", obj_ds, "\n")
    
    objectTADs_posDT <- all_topTADsGenes[[j]][["TADposDT"]]
    objectTADs_posDT$region <- as.character(objectTADs_posDT$region)
    
    obj_g2tDT <- all_topTADsGenes[[j]][["gene2tadDT"]]
    obj_g2tDT$region <- as.character(obj_g2tDT$region)
    # added 30.01: filter the TAD to the ones in g2t -> I want only those with genes inside
    objectTADs_posDT <- objectTADs_posDT[objectTADs_posDT$region %in% obj_g2tDT$region, ]
    stopifnot(nrow(objectTADs_posDT) > 0)
    
    object_allTADs_IR <- IRanges(start = objectTADs_posDT$start, 
                                 width = (objectTADs_posDT$end - objectTADs_posDT$start + 1), 
                                 names=objectTADs_posDT$region)
    
    object_allTADs_GR <- GRanges(ranges = object_allTADs_IR,
                                   seqnames=gsub("(chr.+)_TAD.+", "\\1", objectTADs_posDT$region))
    
    
    TADoverlap_hits <- findOverlaps(query=query_signifTADs_GR, 
                               subject=object_allTADs_GR)
    # CHANGED 30.01 _IR to _GR
    TADoverlaps <- pintersect(query_signifTADs_GR[queryHits(TADoverlap_hits)], 
                           object_allTADs_GR[subjectHits(TADoverlap_hits)])
    
    queryTADs <- names(query_signifTADs_IR[queryHits(TADoverlap_hits)])
    objectTADs <- names(object_allTADs_IR[subjectHits(TADoverlap_hits)])
    stopifnot( length(queryTADs) == length(objectTADs) )
    
    # ensure the matching is done intrachromosomally
    queryTADs_chr <- gsub("(chr.+)_.+", "\\1", queryTADs)
    objectTADs_chr <- gsub("(chr.+)_.+", "\\1", objectTADs)
    stopifnot( queryTADs_chr == objectTADs_chr)
    
    # retrieve the number of common genes for each match
    # take the reference TADs that have a match
    queryTADs_withMatch <- unique(queryTADs)
    
    txt <- paste0("... signif. queryTADs with overlapping objectTADs:\t")
    printAndLog(txt, logFile)
    txt <- paste0(length(queryTADs_withMatch), "/", length(query_signifTADs_GR), "\n")
    printAndLog(txt, logFile)
    
    cat("... retrieve genes for the TADs that have matching object TADs\n")
    
    
    # for each TAD, 1) retrieve the genes that belong to it
    queryTADs_withMatch_genes <- lapply(queryTADs_withMatch, function(curr_tad){
      stopifnot( curr_tad %in% ref_g2tDT$region )
      as.character(ref_g2tDT$entrezID[ref_g2tDT$region == curr_tad])
    })
    names(queryTADs_withMatch_genes) <- queryTADs_withMatch
    
    # queryTADs_withMatch <- queryTADs_withMatch[1:3]
    
    # for each TAD, 2) retrieve the genes of that TAD that matches with it
    
    cat("... retrieve the genes of the matching TADs\n")
    
    queryTADs_withMatch_matchingObjectTADs <- lapply(queryTADs_withMatch, function(curr_tad){
          
          ref_genes <- queryTADs_withMatch_genes[[curr_tad]]
          
          matching_objectTADs <- objectTADs[queryTADs == curr_tad]
          stopifnot(length(matching_objectTADs) > 0)
          
          matching_objectTADs_genes <- lapply(matching_objectTADs, function(obj_tad){
            stopifnot( obj_tad %in% obj_g2tDT$region )
            genes_in_matchingTADs <- as.character(obj_g2tDT$entrezID[obj_g2tDT$region == obj_tad])
            nbr_intersectGenes_in_matchingTADs <- sum(genes_in_matchingTADs %in% ref_genes)
            list(matchingTADs_genes = genes_in_matchingTADs,
                 matchingTADs_nIntersect = nbr_intersectGenes_in_matchingTADs)
          })
          names(matching_objectTADs_genes) <- matching_objectTADs
          matching_objectTADs_genes
        
          # list(matching_objectTADs=matching_objectTADs,
          #      matching_objectTADs_genes=matching_objectTADs_genes)
    })
    names(queryTADs_withMatch_matchingObjectTADs) <- queryTADs_withMatch
    
    ### VERSION 1 -> MATCHING TAD IS THE ONE WITH MOST INTERSECT GENES
    
    # e.g. for the query chr10_TAD1 -> match 2 objectTADs
    # > queryTADs_withMatch_matchingObjectTADs[["chr10_TAD1"]]
    # $chr10_TAD1
    # $chr10_TAD1$matchingTADs_genes
    # [1] "347688"    "439945"    "10771"     "100421369"
    # 
    # $chr10_TAD1$matchingTADs_nIntersect
    # [1] 4
    # 
    # 
    # $chr10_TAD2
    # $chr10_TAD2$matchingTADs_genes
    # [1] "22982"     "100847086" "414235"   
    # 
    # $chr10_TAD2$matchingTADs_nIntersect
    # [1] 3
    
    # > queryTADs_withMatch_matchingObjectTADs[["chr11_TAD130"]]
    # $chr11_TAD108
    # $chr11_TAD108$matchingTADs_genes
    # [1] "246330"    "101928069" "10072"     "582"       "254359"    "89"       
    # [7] "8722"      "55231"     "9973"      "10432"     "100526737" "5936"     
    # [13] "83759"    
    # 
    # $chr11_TAD108$matchingTADs_nIntersect
    # [1] 13
    # 
    # 
    # $chr11_TAD109
    # $chr11_TAD109$matchingTADs_genes
    # [1] "6712"      "79703"     "100422299" "100462788"
    # 
    # $chr11_TAD109$matchingTADs_nIntersect
    # [1] 4
    
    # Filter to retain only the highest intersect
    
    cat("... filter retain highest gene intersect\n")
    
    queryTADs_withMatch_maxInterObjectTADs <- lapply(queryTADs_withMatch_matchingObjectTADs, function(curr_tad){
      tmp_nIntersect <- unlist(lapply(curr_tad, function(x) as.numeric(x[["matchingTADs_nIntersect"]])))
      stopifnot( ! is.na(tmp_nIntersect) )
      obj_maxIntersect_names <- which.max(tmp_nIntersect)
      obj_maxIntersect_names
    })
    stopifnot(names(queryTADs_withMatch_maxInterObjectTADs) == queryTADs_withMatch)
    
    maxIntersectMatchingTAD <- as.character(unlist(
      lapply(queryTADs_withMatch_maxInterObjectTADs, function(x) names(x))))
    
    nIntersectMatchingTAD <- unlist(lapply(queryTADs_withMatch_matchingObjectTADs, function(curr_tad){
      tmp_nIntersect <- unlist(lapply(curr_tad, function(x) as.numeric(x[["matchingTADs_nIntersect"]])))
      stopifnot( ! is.na(tmp_nIntersect) )
      obj_maxIntersect_nbr <- max(tmp_nIntersect)
      obj_maxIntersect_nbr
    }))

    stopifnot(!is.na(nIntersectMatchingTAD))
    stopifnot(is.numeric(nIntersectMatchingTAD))
    
    stopifnot( length(maxIntersectMatchingTAD) == length(queryTADs_withMatch) )
    
    query_nbrGenes <- sapply(queryTADs_withMatch, function(x) 
                          sum(ref_g2tDT$region == x))
    
    object_nbrGenes <- sapply(maxIntersectMatchingTAD, function(x) 
      sum(obj_g2tDT$region == x))
    
    query_matching_ratio <- nIntersectMatchingTAD/query_nbrGenes
    stopifnot(query_matching_ratio >= 0 & query_matching_ratio <= 1 )
    
    maxIntersect_matchDT <- data.frame(
      
      expr_ds = exprds,
      
      query_ds = ref_ds,
      matching_ds = obj_ds,
      
      queryTAD = queryTADs_withMatch,
      query_nbrGenes = query_nbrGenes,
      
      matchingTAD_maxIntersect = maxIntersectMatchingTAD,
      matching_nbrGenes = object_nbrGenes,
        
      query_matching_intersectGenes = nIntersectMatchingTAD,
      query_matching_ratio = query_matching_ratio,
      
      
      stringsAsFactors = FALSE
    )
    stopifnot(maxIntersect_matchDT$query_matching_intersectGenes <= maxIntersect_matchDT$matching_nbrGenes)
    stopifnot(maxIntersect_matchDT$query_matching_intersectGenes <= maxIntersect_matchDT$query_nbrGenes)
    
    maxIntersect_matchDT$queryTAD <- as.character(maxIntersect_matchDT$queryTAD)
    maxIntersect_matchDT <- maxIntersect_matchDT[order(maxIntersect_matchDT$queryTAD),]
    
                # outFile <- file.path(outFolder, paste0(exprds, "_", ref_ds, "_", obj_ds, "_maxIntersect_matchDT.txt"))
                # cat("... writte table with max gene intersect\n")
                # write.table(maxIntersect_matchDT, col.names=TRUE, row.names=FALSE, sep="\t", quote=F, file = outFile)
                # cat(paste0("... written: ", outFile, "\n"))
        
    ### VERSION 2 -> MATCHING TAD IS THE ONE WITH HIGHEST OVERLAP
        
    # TADoverlap_hits <- findOverlaps(query=query_signifTADs_GR, 
    #                                 subject=object_allTADs_GR)
    # 
    # TADoverlaps <- pintersect(query_signifTADs_GR[queryHits(TADoverlap_hits)], 
    #                           object_allTADs_GR[subjectHits(TADoverlap_hits)])
    percentTADoverlap <- width(TADoverlaps)/width(query_signifTADs_GR[queryHits(TADoverlap_hits)])
    stopifnot( percentTADoverlap > 0 & percentTADoverlap <= 1 )
    

    overlap_matchDT <- data.frame(
      
      expr_ds = exprds,
      
      query_ds = ref_ds,
      matching_ds = obj_ds,
      
      queryTAD = names(query_signifTADs_GR[queryHits(TADoverlap_hits)]),
      query_size = width(query_signifTADs_GR[queryHits(TADoverlap_hits)]),
      
      matchingTAD = names(object_allTADs_GR[subjectHits(TADoverlap_hits)]),
      matching_size = width(object_allTADs_GR[subjectHits(TADoverlap_hits)]),
      
      query_matching_overlap = percentTADoverlap,
      stringsAsFactors = FALSE
    )
    


        
    maxOverlap_matchDT <- do.call(rbind, by(data = overlap_matchDT, INDICES = overlap_matchDT$queryTAD, FUN=function(x) {
      x[which.max(x$query_matching_overlap),]
    }))

    colnames(maxOverlap_matchDT)[  colnames(maxOverlap_matchDT) == "matchingTAD"] <- "matchingTAD_maxOverlap"
    
    maxOverlap_matchDT$queryTAD <- as.character(maxOverlap_matchDT$queryTAD)
    maxOverlap_matchDT <- maxOverlap_matchDT[order(maxOverlap_matchDT$queryTAD),]
    
                # outFile <- file.path(outFolder, paste0(exprds, "_", ref_ds, "_", obj_ds, "_maxOverlap_matchDT.txt"))
                # cat("... writte table with max overlap\n")
                # write.table(maxOverlap_matchDT, col.names=TRUE, row.names=FALSE, sep="\t", quote=F, file = outFile)
                # cat(paste0("... written: ", outFile, "\n"))
    
    
    # list(
    #   maxIntersect = maxIntersect_matchDT,
    #   maxOverlap = maxOverlap_matchDT
    # )
    #break
    
    # expr_ds	query_ds	matching_ds	queryTAD	query_nbrGenes	matchingTAD_maxIntersect	matching_nbrGenes	query_matching_intersectGenes
    # expr_ds	query_ds	matching_ds	queryTAD	query_size	matchingTAD_maxOverlap	matching_size	query_matching_overlap
    
    maxIntersectOverlap_matchDT <- merge(maxIntersect_matchDT, maxOverlap_matchDT, by=c("expr_ds", "query_ds", "matching_ds", "queryTAD"))
    maxIntersectOverlap_matchDT
    
  } # end iterating over the j [object dataset] -> ref_maxIntersectOverlap_matchDT
  #break
  ref_maxIntersectOverlap_matchDT
} # end iterating over the i [reference dataset] -> all_maxIntersectOverlap_matchDT

outFile <- file.path(outFolder, paste0(exprds, "_", ref_ds, "_", obj_ds, "_top", topThresh, "_all_maxIntersectOverlap_matchDT.Rdata"))
cat("... writte table with max overlap and intersect\n")
save(all_maxIntersectOverlap_matchDT,file = outFile)
cat(paste0("... written: ", outFile, "\n"))


outDT <- all_maxIntersectOverlap_matchDT

# outDT <- outDT[order(outDT[, "query_matching_ratio"]),]


outDT[, "query_matching_ratio"] <- round(outDT[, "query_matching_ratio"], 4)
outDT[, "query_matching_overlap"] <- round(outDT[,"query_matching_overlap"], 4)


outFile <- file.path(outFolder, paste0(exprds, "_", ref_ds, "_", obj_ds, "_top", topThresh, "_all_maxIntersectOverlap_matchDT.txt"))
cat("... writte table with max overlap and intersect\n")
write.table(outDT, col.names=TRUE, row.names=FALSE, sep="\t", quote=F, file = outFile)
cat(paste0("... written: ", outFile, "\n"))

shortDT <- outDT 
shortDT <- shortDT[, c("expr_ds", "query_ds", "matching_ds", "queryTAD", "matchingTAD_maxIntersect", "query_matching_ratio")]
shortDT <- shortDT[order(shortDT[, "query_matching_ratio"], decreasing = TRUE),]
outFile <- file.path(outFolder, paste0(exprds, "_", ref_ds, "_", obj_ds, "_top", topThresh, "_all_maxIntersectOverlap_matchDT_shortDT.txt"))
cat("... writte table with max overlap and intersect\n")
write.table(shortDT, col.names=TRUE, row.names=FALSE, sep="\t", quote=F, file = outFile)
cat(paste0("... written: ", outFile, "\n"))


cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))



# ######################################################################################
# cat(paste0("... written: ", outFile, "\n"))
# cat(paste0("... written: ", listFile, "\n"))





