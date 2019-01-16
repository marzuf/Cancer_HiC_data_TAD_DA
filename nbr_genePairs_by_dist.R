
hicds <- "ENCSR079VIJ_G401_40kb"
exprds <- "TCGAkich_norm_kich"

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2)

stopifnot(dir.exists(hicds))

# retrieve the information from the loess model, so no need to recompute everything !
# (loess model computed in:  Rscript AUC_coexprDist_withFam_sortNoDup_otherTADfile_otherFamFile.R)

geneFile <- file.path("PIPELINE/OUTPUT_FOLDER", hicds, exprds, "0_prepGeneData/pipeline_geneList.Rdata")
stopifnot(file.exists(geneFile))
geneList <- eval(parse(text = load(geneFile)))
# in geneList, entrezID is the vector items (not the names)
geneList <- as.character(geneList)
# ENSG00000000003 ENSG00000000005 ENSG00000000419 ENSG00000000457 ENSG00000000460 
# "7105"         "64102"          "8813"         "57147"         "55732" 

distFile <- file.path("CREATE_DIST_SORTNODUP", hicds, "all_dist_pairs.Rdata")
stopifnot(file.exists(distFile))
distDT <- eval(parse(text = load(distFile)))
distDT$gene1 <- as.character(distDT$gene1)
distDT$gene2 <- as.character(distDT$gene2)

# take only the genes for the given dataset
distDT <- distDT[distDT$gene1 %in% geneList & distDT$gene2 %in% geneList,] 
