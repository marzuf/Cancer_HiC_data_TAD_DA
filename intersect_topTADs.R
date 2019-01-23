

# 
# Rscript intersect_topTADs.R TCGAluad_mutKRAS_mutEGFR 10 ENCSR444WCZ_A549_40kb NCI-H460_40kb pipelineConsensus
# Rscript intersect_topTADs.R TCGAluad_mutKRAS_mutEGFR 10 ENCSR444WCZ_A549_40kb
# Rscript intersect_topTADs.R TCGAluad_mutKRAS_mutEGFR 10 NCI-H460_40kb
# Rscript intersect_topTADs.R TCGAluad_mutKRAS_mutEGFR 10 pipelineConsensus
# 
# Rscript intersect_topTADs.R TCGAkich_norm_kich 10 ENCSR079VIJ_G401_40kb ENCSR401TBQ_Caki2_40kb ENCSR079VIJ_G401ENCSR401TBQ_Caki2_40kb pipelineConsensus
# Rscript intersect_topTADs.R TCGAkich_norm_kich 10 ENCSR079VIJ_G401_40kb
# Rscript intersect_topTADs.R TCGAkich_norm_kich 10 ENCSR401TBQ_Caki2_40kb
# Rscript intersect_topTADs.R TCGAkich_norm_kich 10 ENCSR079VIJ_G401ENCSR401TBQ_Caki2_40kb
# Rscript intersect_topTADs.R TCGAkich_norm_kich 10 pipelineConsensus
# 
# Rscript intersect_topTADs.R TCGAskcm_wt_mutCTNNB1 10 ENCSR312KHQ_SK-MEL-5_40kb ENCSR862OGI_RPMI-7951_40kb ENCSR312KHQ_SK-MEL-5ENCSR862OGI_RPMI-7951_40kb pipelineConsensus
# Rscript intersect_topTADs.R TCGAskcm_wt_mutCTNNB1 10 ENCSR312KHQ_SK-MEL-5_40kb
# Rscript intersect_topTADs.R TCGAskcm_wt_mutCTNNB1 10 ENCSR862OGI_RPMI-7951_40kb
# Rscript intersect_topTADs.R TCGAskcm_wt_mutCTNNB1 10 ENCSR312KHQ_SK-MEL-5ENCSR862OGI_RPMI-7951_40kb
# Rscript intersect_topTADs.R TCGAskcm_wt_mutCTNNB1 10 pipelineConsensus
# 

startTime <- Sys.time()

require(foreach)

setDir=""
setDir <- "~/media/electron"
setDir="" 

outFolder <- file.path("INTERSECT_topTADs")
dir.create(outFolder, recursive = TRUE)

args <- commandArgs(trailingOnly = TRUE)

all_hicds <- c("ENCSR312KHQ_SK-MEL-5_40kb", "ENCSR862OGI_RPMI-7951_40kb", "ENCSR312KHQ_SK-MEL-5ENCSR862OGI_RPMI-7951_40kb")
exprds <- "TCGAskcm_wt_mutCTNNB1"
topThresh <- 5
stopifnot(length(args) > 2)


all_hicds <- args[3:length(args)]
exprds <- args[1]
nTop <- as.numeric(args[2])

stopifnot(!is.na(nTop))
stopifnot(nTop > 0)

entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
entrezDT <- read.delim(entrezDT_file, header=TRUE, stringsAsFactors = FALSE)
entrezDT$entrezID <- as.character(entrezDT$entrezID)
head(entrezDT)

#PIPELINE/OUTPUT_FOLDER/GSE105318_DLD1_40kb/TCGAcoad_msi_mss//emp_pval_combined.Rdata
script11_name <- "11_runEmpPvalCombined"

all_topTADsGenes <- foreach(hicds = all_hicds) %dopar% {
  
  if(hicds == "pipelineConsensus") {
    pipOutDir <- file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom", "OUTPUT_FOLDER",  exprds)
    
    # file with assignment from entrez to all regions
    gene2tadDT_file <- paste0(setDir, 
                              "/mnt/ed4/marie/gene_data_final/consensus_TopDom_covThresh_r0.6_t80000_v0_w-1_final/genes2tad/all_genes_positions.txt") 
  } else {
    stopifnot(dir.exists(hicds))
    pipOutDir <- file.path("PIPELINE", "OUTPUT_FOLDER", hicds, exprds)
    gene2tadDT_file <- file.path(hicds, "genes2tad/all_genes_positions.txt")
  }
  stopifnot(file.exists(gene2tadDT_file))
  g2tDT <- read.delim(gene2tadDT_file, stringsAsFactors = FALSE, header=F, col.names=c("entrezID", "chromo", "start", "end", "region"))
  
  stopifnot(g2tDT$entrezID %in% entrezDT$entrezID)
  g2t2sDT <- merge(g2tDT[, c("entrezID", "region")], entrezDT[,c("entrezID", "symbol")], by="entrezID", all.x = TRUE, all.y = FALSE )
  head(g2t2sDT)
  stopifnot(!is.na(g2t2sDT))
  g2t2sDT$entrezID <- as.character(g2t2sDT$entrezID)
  g2t2sDT$region <- as.character(g2t2sDT$region)
  g2t2sDT$symbol <- as.character(g2t2sDT$symbol)
  
  stopifnot(dir.exists(pipOutDir))
  stopifnot(dir.exists(file.path(pipOutDir, script11_name)))
  
  tad_pvalFile <- file.path(pipOutDir, script11_name, "emp_pval_combined.Rdata")
  stopifnot(file.exists(tad_pvalFile))
  tad_pvals <- eval(parse(text = load(tad_pvalFile)))
  
  adj_tad_pvals <- sort(p.adjust(tad_pvals, method="BH"))
  
  if(topThresh > 1) {
    topThresh <- min(c(topThresh, length(adj_tad_pvals)))
    pvalThresh <- as.numeric(adj_tad_pvals[topThresh])
    stopifnot(!is.na(pvalThresh))
    stopifnot(pvalThresh <= 1)
  } else {
    pvalThresh <- topThresh
    topThresh <- min(which(adj_tad_pvals <= pvalThresh))
  }
  
  top_pvals <- adj_tad_pvals[adj_tad_pvals <= topThresh]
  stopifnot(!is.na(top_pvals))
  
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
  
  topTADs_genesDT <- merge(topTADsDT, topGenesDT[, c("top_tads", "symbol")], by = "top_tads", all.x=TRUE, all.y=TRUE)
  
  topTADs_genesDT$top_tads <- as.character(topTADs_genesDT$top_tads)
  topTADs_genesDT$symbol <- as.character(topTADs_genesDT$symbol)
  topTADs_genesDT <- topTADs_genesDT[order(topTADs_genesDT$top_pvals, topTADs_genesDT$top_tads, topTADs_genesDT$symbol),]
  head(topTADs_genesDT)
  
  topTADs_genesDT
}
names(all_topTADsGenes) <- all_hicds
str(all_topTADsGenes)

intersectGenes <- Reduce(intersect, lapply(all_topTADsGenes, function(x) x[["symbol"]]))

filtered_topTADsGenes <- lapply(all_topTADsGenes, function(x) x[x[["symbol"]] %in% intersectGenes,, drop=F])
filtered_topTADsGenes <- lapply(filtered_topTADsGenes, function(x) x[order(x[["symbol"]]),])

filtered_topTADsGenesDT <- do.call('cbind', filtered_topTADsGenes)

all_symbolsDT <- filtered_topTADsGenesDT[, grepl("\\.symbol$", colnames(filtered_topTADsGenesDT)), drop=F]
stopifnot(apply(all_symbolsDT, 1, function(x) length(unique(x)) == 1))
all_symbols <- as.character(all_symbolsDT[,1])
filtered_topTADsGenesDT[, grepl("\\.symbol$", colnames(filtered_topTADsGenesDT))] <- NULL
rownames(filtered_topTADsGenesDT) <- all_symbols
head(filtered_topTADsGenesDT)

all_ranksDT <- filtered_topTADsGenesDT[, grepl("\\.tad_rank$", colnames(filtered_topTADsGenesDT)), drop=F]
filtered_topTADsGenesDT$rank_mean <- rowMeans(all_ranksDT)
head(filtered_topTADsGenesDT)
filtered_topTADsGenesDT$gene <- as.character(rownames(filtered_topTADsGenesDT))

filtered_topTADsGenesDT[, grepl("\\.tad_rank$", colnames(filtered_topTADsGenesDT))] <- NULL


filtered_topTADsGenesDT <- filtered_topTADsGenesDT[, c(which(colnames(filtered_topTADsGenesDT) == "gene"),
                                                       which(colnames(filtered_topTADsGenesDT) == "rank_mean"),
                                                       which(! colnames(filtered_topTADsGenesDT) %in% c("gene", "rank_mean"))), drop=F]

if(length(hicds) == 1) {
  filtered_topTADsGenesDT <- filtered_topTADsGenesDT[order(filtered_topTADsGenesDT$rank_mean, filtered_topTADsGenesDT$gene),, drop=F ]
} else {
  filtered_topTADsGenesDT <- filtered_topTADsGenesDT[order(filtered_topTADsGenesDT$gene, filtered_topTADsGenesDT$rank_mean),, drop=F ]
}
head(filtered_topTADsGenesDT)

cat("\n\n\n")

write.table(head(filtered_topTADsGenesDT), col.names=T, row.names=F, sep="\t", quote=F, file="")

outFile <- file.path(outFolder, 
                     paste0(exprds, "_", paste0(hicds, collapse=""), "_", nTop, ".txt"))

write.table(filtered_topTADsGenesDT, col.names=T, row.names=F, sep="\t", quote=F, file=outFile)

cat(paste0("... written: ", outFile, "\n"))
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))

# resultDT <- data.frame(
#   # top_pvals = round(top_pvals, 4),
#   # top_pvals = top_pvals,
#   top_pvals = sprintf("%.2e",top_pvals),
#   top_tads = top_tads,
#   stringsAsFactors = FALSE
# )
# 
# cat(paste0("********** ", hicds, " - ", exprds, "\n\n"))
# cat(paste0("> topThresh\t=\t", topThresh, "\n"))
# cat(paste0("> pvalThresh\t=\t", sprintf("%.2e",pvalThresh), "\n\n"))
# cat(paste0("*** Top ranking TADs (adj. emp. p-val. combined): \n"))
# 
# #write.table(head(resultDT), col.names=TRUE, row.names=FALSE, sep="\t", quote=F, file = "")
# 
# write.table(resultDT, col.names=TRUE, row.names=FALSE, sep="\t", quote=F, file = "")









