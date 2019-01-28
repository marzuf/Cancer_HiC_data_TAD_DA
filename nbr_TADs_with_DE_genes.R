# Rscript nbr_TADs_with_DE_genes.R

cat("> START nbr_TADs_with_DE_genes.R\n")

startTime <- Sys.time()

hicds <- "ENCSR079VIJ_G401_40kb"
exprds <- "TCGAkich_norm_kich"

args <- commandArgs(trailingOnly = TRUE)
hicds <- args[1]
exprds<- args[2]

signifThresh <- 0.05

plotCex <- 1.2

plotType <- "svg"
myHeight <- ifelse(plotType == "png", 400, 7)
myWidth <- ifelse(plotType == "png", 600, 10)

source("utils_fct.R")

outFold <- "NBR_TADS_WITH_DE_GENES"
dir.create(outFold, recursive = TRUE)

SSHFS <- TRUE
if(SSHFS) logFile <- ""

pipFold <- file.path("PIPELINE", "OUTPUT_FOLDER", hicds, exprds)
stopifnot(dir.exists(pipFold))

gene2tad_DT_file <- file.path(hicds, "genes2tad", "all_genes_positions.txt")
stopifnot(file.exists(gene2tad_DT_file))

g2t_DT <- read.delim(gene2tad_DT_file, header=F, stringsAsFactors = FALSE,
                     col.names=c("entrezID", "chromo", "start", "end", "region"))
g2t_DT$entrezID <- as.character(g2t_DT$entrezID) 
head(g2t_DT)

topTable_file <- file.path(pipFold, "1_runGeneDE", "DE_topTable.Rdata")
stopifnot(file.exists(topTable_file))
topTable_DT <- eval(parse(text = load(topTable_file)))
head(topTable_DT)
topTable_DT$genes <- as.character(topTable_DT$genes)
# topTable_DT$genes might be ENSEMBL, entrez, symbol, etc.
# geneList: values are entrez, names can be ENSEMBL, entrez, etc.

geneList_file <- file.path(pipFold, "0_prepGeneData", "pipeline_geneList.Rdata")
stopifnot(file.exists(geneList_file))
geneList <- eval(parse(text = load(geneList_file)))

# should be tested in this direction because DE analysis was done with all available genes !
stopifnot(names(geneList) %in% topTable_DT$genes)

stopifnot(geneList %in% g2t_DT$entrezID)

topTable_DT <- topTable_DT[topTable_DT$genes %in% names(geneList),]
topTable_DT$entrezID <- sapply(topTable_DT$genes, function(x) {
  eID <- geneList[as.character(x)]
  stopifnot(length(eID) == 1)
  stopifnot(!is.na(eID))
  as.character(eID)
})
stopifnot(topTable_DT$entrezID %in% g2t_DT$entrezID)

head(topTable_DT)
head(g2t_DT)

DE_g2t_DT <- merge(topTable_DT[, c("entrezID", "adj.P.Val")], g2t_DT[, c("entrezID", "region")], by="entrezID")

nbr_tads_tot <- length(unique(DE_g2t_DT$region))
txt <- paste0("... total nbr of TADs:\t", nbr_tads_tot, "\n")
printAndLog(txt, logFile)

# observed number of TADs containing DE genes
obs_nbr_tads_withDE <- length(unique(DE_g2t_DT$region[DE_g2t_DT$adj.P.Val <= signifThresh]))
txt <- paste0("... nbr of TADs with DE genes:\t", obs_nbr_tads_withDE, "\n")
printAndLog(txt, logFile)


# compute the same values for the permutation data
permutDT_file <- file.path(pipFold, "5_runPermutationsMedian", "permutationsDT.Rdata")
stopifnot(file.exists(permutDT_file))
cat("... load permutDT \n")
permutDT <- eval(parse(text = load(permutDT_file)))

DE_genes <- DE_g2t_DT$entrezID[DE_g2t_DT$adj.P.Val <= signifThresh]
stopifnot(DE_genes %in% rownames(permutDT))

permutDT <- permutDT[rownames(permutDT) %in% DE_genes,]

perm_nbr_tads_withDE <- apply(permutDT, 2, function(x) {length(unique(x))})


save(perm_nbr_tads_withDE, file = file.path(outFold, "perm_nbr_tads_withDE.Rdata"))

stopifnot(!is.na(perm_nbr_tads_withDE))

outFile <- file.path(outFold, paste0(hicds, "_", exprds, "_TADs_with_DE_genes.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot(density(perm_nbr_tads_withDE),
     main = paste0(hicds, " - ", exprds, " - # TADs with DE genes (adj. p-val. <= ", signifThresh, ")"), 
     xlab = "# TADs with DE genes",
     cex.axis=plotCex, cex.lab=plotCex)
abline(v = obs_nbr_tads_withDE, col = "red", lty = 2, lwd=1.5 )
legend("topright",
       c("obs. value", "dist. perm. values"),
       lty = c(2,1),
       col = c("red", "black"),
       bty="n")
mtext(side=3, paste0("dist. permut. (n=", length(perm_nbr_tads_withDE), ") vs. obs. value (", obs_nbr_tads_withDE, "/", nbr_tads_tot,")"))
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))
###################################################################
cat(paste0("***DONE\n", startTime,"\n", Sys.time(), "\n"))    



