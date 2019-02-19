
# Rscript intersect_topTADs_acrossDS.R 3

startTime <- Sys.time()

cat("> START intersect_topTADs_acrossDS.R \n")

SSHFS <- FALSE

require(foreach)
require(doMC)
require(IRanges)
require(GenomicRanges)
source("utils_fct.R")

printVar <- function(x){
  cat(paste0(x, " = ", eval(parse(text=x)), "\n"))
}

registerDoMC(ifelse(SSHFS, 2, 40))

build_signifTADs_allDS_data <- TRUE


setDir <- ifelse(SSHFS, "~/media/electron", "")

topThresh <- 3

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) > 0)
topThresh <- as.numeric(args[1])

if(topThresh == 1) {
  warning("topThresh == 1 is ambiguous; will be considered as a nTop not pval thresh !\n")
}

outFolder <- file.path("INTERSECT_topTADs_ACROSSDS", paste0("top", topThresh))
dir.create(outFolder, recursive = TRUE)


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

pipOutFolder <- file.path("PIPELINE", "OUTPUT_FOLDER")

all_hicexpr_ds <- unname(unlist(sapply(list.files(pipOutFolder, full.names = TRUE), function(x) file.path(basename(x),  list.files(x)))))
stopifnot(dir.exists(file.path(pipOutFolder, all_hicexpr_ds)))


ds=all_hicexpr_ds[1]

if(build_signifTADs_allDS_data){
  
cat("... start building signifTADs_allDS_data \n")
signifTADs_allDS_data <- foreach(ds = all_hicexpr_ds) %dopar% {
  
  cat("... start: ", ds, "\n")
  
  hicds <- dirname(ds)
  exprds <- basename(ds)
  stopifnot(dir.exists(hicds))
  dsPipOutDir <- file.path(pipOutFolder, ds)
  stopifnot(dir.exists(dsPipOutDir))
  gene2tadDT_file <- file.path(hicds, "genes2tad/all_genes_positions.txt")
  TADposDT_file <- file.path(hicds, "genes2tad/all_assigned_regions.txt")
  stopifnot(file.exists(gene2tadDT_file))
  stopifnot(file.exists(TADposDT_file))
  
  # RETRIEVE SIGNIF. TADs 
  stopifnot(dir.exists(file.path(dsPipOutDir, script11_name)))
  tad_pvalFile <- file.path(dsPipOutDir, script11_name, "emp_pval_combined.Rdata")
  stopifnot(file.exists(tad_pvalFile))
  tad_pvals <- eval(parse(text = load(tad_pvalFile)))
  adj_tad_pvals <- sort(p.adjust(tad_pvals, method="BH"))
  
  # RETRIEVE PIPELINE GENES
  # subset the genes used in the pipeline
  stopifnot(dir.exists(file.path(dsPipOutDir, script0_name)))
  # geneList: values are entrez, names can be ENSEMBL, entrez, etc.
  geneList_file <- file.path(dsPipOutDir, script0_name, "pipeline_geneList.Rdata")
  stopifnot(file.exists(geneList_file))
  geneList <- eval(parse(text = load(geneList_file)))
  
  
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
  
  # RETRIEVE RANKS OF SIGNIF TADs (-> signifDT)
  tad_ranks <- rank(top_pvals, ties="min")
  stopifnot(top_tads %in% names(tad_ranks))
  tad_ranks <- tad_ranks[top_tads]
  
  # RETRIEVE POSITION OF SIGNIF TADs (-> posDT)
  TADposDT <- read.delim(TADposDT_file, stringsAsFactors = FALSE, header=F, col.names=c("chromo","region", "start", "end"))
  stopifnot(is.numeric(TADposDT$start))
  TADposDT <- TADposDT[TADposDT$region %in% top_tads,]
  TADposDT <- TADposDT[match(top_tads, TADposDT$region),]
  stopifnot(nrow(TADposDT) > 0)
  stopifnot(!is.na(TADposDT))
  stopifnot(TADposDT$region %in% top_tads)
  stopifnot(top_tads %in% TADposDT$region)
  
  stopifnot(TADposDT$region == top_tads)
  head(TADposDT)
  
  
  # RETRIEVE GENES ID AND GENE SYMBOLS FROM SIGNIF TADs (-> geneDT)
  
  
  g2tDT <- read.delim(gene2tadDT_file, stringsAsFactors = FALSE, header=F, col.names=c("entrezID", "chromo", "start", "end", "region"))
  stopifnot(geneList %in% g2tDT$entrezID)
  stopifnot(top_tads %in% g2tDT$region)
  g2tDT <- g2tDT[g2tDT$entrezID %in% geneList &
                   g2tDT$region %in% top_tads,]
  stopifnot(g2tDT$entrezID %in% entrezDT$entrezID)
  g2t2sDT <- merge(g2tDT[, c("entrezID", "region")], entrezDT[,c("entrezID", "symbol")], by="entrezID", all.x = TRUE, all.y = FALSE )
  stopifnot(!is.na(g2t2sDT))
  g2t2sDT$entrezID <- as.character(g2t2sDT$entrezID)
  g2t2sDT$region <- as.character(g2t2sDT$region)
  g2t2sDT$symbol <- as.character(g2t2sDT$symbol)
  head(g2t2sDT)
  stopifnot(!duplicated(g2t2sDT$entrezID))
  stopifnot(top_tads %in% g2t2sDT$region)
  head(g2t2sDT)
  #geneDT_tmp <- g2t2sDT[order(match(top_tads, g2t2sDT$region), g2t2sDT$entrezID),] # will not work because of the  duplicated values in g2t2sDT$region !!!
  
  geneDT_tmp <- g2t2sDT[order(unlist(sapply(g2t2sDT$region, function(x) which(top_tads == x) )), as.character(g2t2sDT$entrezID)),]
  
  stopifnot(nrow(geneDT_tmp) > 0)
  stopifnot(!is.na(geneDT_tmp))
  stopifnot(geneDT_tmp$region %in% top_tads)
  head(geneDT_tmp)
  stopifnot(top_tads  %in% geneDT_tmp$region)

  geneDT <- data.frame(
    ID = paste(hicds, exprds, as.character(geneDT_tmp$region), sep="_"),
    GENE = as.character(geneDT_tmp$entrezID), 
    SYMBOL = as.character(geneDT_tmp$symbol),
    stringsAsFactors = FALSE
    )
  
  
  
  # RETRIEVE MATCHING TADs
  
  # BUIlD THE TABLES:
  # (matchDT built later)
  id_col <- paste(hicds, exprds, top_tads, sep="_")
  
  idDT <- data.frame(
    ID = id_col,
    HICDS = hicds,
    EXPRDS = exprds,
    TAD = top_tads,
    stringsAsFactors = FALSE
  )  
  rownames(idDT) <- NULL
  head(idDT)
  
  signifDT <- data.frame(
    ID = id_col,
    RANK = tad_ranks,              # not sure necessary: could be retrieve rank(PVAL)
    PVAL = as.numeric(top_pvals),
    stringsAsFactors = FALSE
  )
  rownames(signifDT) <- NULL
  head(signifDT)
  
  stopifnot(TADposDT$region == top_tads)
  posDT <- data.frame(
    ID = id_col,
    CHROMO = TADposDT$chromo,
    START  = TADposDT$start,
    END = TADposDT$end,
    stringsAsFactors = FALSE
  )
  rownames(posDT) <- NULL
  head(posDT)
  
  stopifnot(id_col == unique(geneDT$ID))
  
  matchDT <- data.frame(
    ID = id_col,
    MATCHID = NA_character_,
    MATCHRATIO = NA_real_,
    stringsAsFactors = FALSE
  )
  
  list(
    ids=id_col,
    idDT = idDT,
    geneDT = geneDT,
    signifDT=signifDT,
    posDT=posDT,
    matchDT=matchDT
  )
}
names(signifTADs_allDS_data) <- all_hicexpr_ds
outFile <- file.path(outFolder, "signifTADs_allDS_data.Rdata")
save(signifTADs_allDS_data, file = outFile)
cat(paste0("... written: ", outFile, "\n"))


} else { # end-if(build_signifTADs_allDS_data)
  outFile <- file.path(outFolder, "signifTADs_allDS_data.Rdata")
  signifTADs_allDS_data <- eval(parse(text = load(outFile)))
}

load("INTERSECT_topTADs_ACROSSDS/signifTADs_allDS_data.Rdata")
head(signifTADs_allDS_data[["NCI-H460_40kb/TCGAlusc_norm_lusc"]][["geneDT"]])
head(signifTADs_allDS_data[["NCI-H460_40kb/TCGAlusc_norm_lusc"]][["signifDT"]])
head(signifTADs_allDS_data[["NCI-H460_40kb/TCGAlusc_norm_lusc"]][["posDT"]])
head(signifTADs_allDS_data[["NCI-H460_40kb/TCGAlusc_norm_lusc"]][["matchDT"]])
head(signifTADs_allDS_data[["NCI-H460_40kb/TCGAlusc_norm_lusc"]][["idDT"]])


all_geneDT <- foreach(ds = names(signifTADs_allDS_data), .combine='rbind') %dopar% {
  signifTADs_allDS_data[[paste0(ds)]][["geneDT"]]
}
printVar("nrow(all_geneDT)")

all_signifDT <- foreach(ds = names(signifTADs_allDS_data), .combine='rbind') %dopar% {
  signifTADs_allDS_data[[paste0(ds)]][["signifDT"]]
}
printVar("nrow(all_signifDT)")


# now I have to build matchDT ...
ds=names(signifTADs_allDS_data)[1]
inputGR_DT <- foreach(ds = names(signifTADs_allDS_data), .combine='rbind') %dopar% {
  posDT <- signifTADs_allDS_data[[paste0(ds)]][["posDT"]]
  idDT <- signifTADs_allDS_data[[paste0(ds)]][["idDT"]]
  outDT <- merge(x=posDT, y=idDT, by="ID", all.x=TRUE, all.y=TRUE )
  head(outDT)
  stopifnot(!is.na(outDT))
  outDT
}

inputGR_DT$ID <- as.character(inputGR_DT$ID)
inputGR_DT$CHROMO <- as.character(inputGR_DT$CHROMO)
inputGR_DT$HICDS <- as.character(inputGR_DT$HICDS)
inputGR_DT$EXPRDS <- as.character(inputGR_DT$EXPRDS)
stopifnot(is.numeric(inputGR_DT$START))
stopifnot(is.numeric(inputGR_DT$END))
head(inputGR_DT)

nAllTADs <- nrow(inputGR_DT)
printVar("nrow(inputGR_DT)")

query_IR <- IRanges(start = inputGR_DT$START, 
                  width = (inputGR_DT$END - inputGR_DT$START + 1), 
                  names=inputGR_DT$ID)

query_GR <- GRanges(ranges = query_IR,
                    seqnames=inputGR_DT$CHROMO)
mcols(query_GR)$ID <- inputGR_DT$ID
stopifnot( names(query_GR) == mcols(query_GR)$ID )
head(mcols(query_GR))

cat("... compute overlap among all TADs")
IDoverlap_hits_all <- findOverlaps(query=query_GR,
                           subject=query_GR)
IDoverlap_hits <- IDoverlap_hits_all[queryHits(IDoverlap_hits_all) != subjectHits(IDoverlap_hits_all)] 

# remove matching with itself
printVar("length(IDoverlap_hits_all)")
printVar("length(IDoverlap_hits)")
stopifnot( (length(IDoverlap_hits_all)-length(IDoverlap_hits) ) == length(query_IR) )

# !!! pintersect removes the match with itself !
# pintersect(query_IR["ENCSR346DCU_LNCaP_40kb_TCGAprad_norm_prad_chr1_TAD258",
#                     "ENCSR346DCU_LNCaP_40kb_TCGAprad_norm_prad_chr1_TAD258"], 
#            query_IR["ENCSR346DCU_LNCaP_40kb_TCGAprad_norm_prad_chr1_TAD258",
#                     "ENCSR079VIJ_G401_40kb_TCGAkich_norm_kich_chr1_TAD274"])
queryID <- names(query_GR[queryHits(IDoverlap_hits)])
matchingID <- names(query_GR[subjectHits(IDoverlap_hits)])
IDsOverlapDT <- data.frame(
  queryID = queryID,
  matchingID = matchingID,
  overlapBP = width(pintersect(query_IR[queryID], query_IR[matchingID])),
    stringsAsFactors = FALSE)
head(IDsOverlapDT)

# check that the matching was done intrachromosomally only !
queryChromo <- gsub(".+_(.+?)_TAD.+?$","\\1", IDsOverlapDT$queryID)
matchingChromo <- gsub(".+_(.+?)_TAD.+?$","\\1", IDsOverlapDT$matchingID)
stopifnot(queryChromo == matchingChromo)

printVar("nAllTADs")
all_query_IDs <- unique(IDsOverlapDT$queryID)
printVar("length(all_query_IDs)")

query_id <- all_query_IDs[1]

head(all_geneDT)
stopifnot(is.character(all_geneDT$ID))
stopifnot(is.character(all_geneDT$GENE))
stopifnot(is.character(all_geneDT$SYMBOL))

cat("... start iterating over all queryIDs to build all_matchDT_noNA\n")

all_matchDT_noNA <- foreach(query_id = all_query_IDs, .combine='rbind') %dopar% {
  
  cat("...... start ", query_id, "\n")
  
  query_genes <- all_geneDT$GENE[all_geneDT$ID == query_id]
  stopifnot(length(query_genes) > 0)
  
  all_matching_IDs <- IDsOverlapDT$matchingID[IDsOverlapDT$queryID == query_id]
  nMatches <- length(all_matching_IDs)
  stopifnot(nMatches > 0)
  
  all_matching_ratio_unstd <- sapply(all_matching_IDs, function(matching_id) {
    matching_genes <- all_geneDT$GENE[all_geneDT$ID == matching_id]
    stopifnot(length(matching_genes) > 0)
    intersectGenes <- intersect(matching_genes, query_genes)
    unionGenes <- union(matching_genes, query_genes)
    length(intersectGenes)/length(unionGenes)
   })
  stopifnot(names(all_matching_ratio_unstd) == all_matching_IDs)
  all_matching_ratio <- sort(all_matching_ratio_unstd, decreasing = TRUE)
  stopifnot( length(all_matching_ratio) == nMatches)
  
  query_hicds <- inputGR_DT$HICDS[inputGR_DT$ID == query_id]
  query_exprds <- inputGR_DT$EXPRDS[inputGR_DT$ID == query_id]
  query_tad <- inputGR_DT$TAD[inputGR_DT$ID == query_id]
  stopifnot( paste0(query_hicds, "_", query_exprds, "_", query_tad) == query_id)
  
  sortedIDs <- names(all_matching_ratio)
  stopifnot(!duplicated(sortedIDs))
  
  matching_hicds <- as.character(sapply(sortedIDs, function(matching_id) {
    inputGR_DT$HICDS[inputGR_DT$ID == matching_id]
  }))
  stopifnot( length(matching_hicds) == length(sortedIDs) )
  
  matching_exprds <- as.character(sapply(sortedIDs, function(matching_id) {
    inputGR_DT$EXPRDS[inputGR_DT$ID == matching_id]
  }))
  stopifnot( length(matching_exprds) == length(sortedIDs) )

   matching_tads <- as.character(sapply(sortedIDs, function(matching_id) {
     inputGR_DT$TAD[inputGR_DT$ID == matching_id]
   }))
   stopifnot( length(matching_tads) == length(sortedIDs) )
   
  all_matching_ratio <- as.numeric(all_matching_ratio)
  stopifnot(!is.na(all_matching_ratio))
 
  matchDT <- data.frame(
    query_id = query_id,
   query_hicds = rep(query_hicds, nMatches),
   query_exprds = rep(query_exprds, nMatches),
   queryTAD = rep(query_tad, nMatches),
   
   matching_id = sortedIDs,
   matching_hicds = matching_hicds,
   matching_exprds = matching_exprds,
   matchingTAD = matching_tads,
   
   matchingRatio = all_matching_ratio,
   
   stringsAsFactors = FALSE
  )
}

# ADD NA FOR THE ONES WITHOUT ANY MATCH
all_ds_tads <- as.character(inputGR_DT$ID)
stopifnot( length(all_ds_tads) == nAllTADs)

noMatch_tads <- all_ds_tads[!all_ds_tads %in% all_matchDT_noNA$query_id]
printVar("nAllTADs")
printVar("length(noMatch_tads)")

cat("... start iterating over IDs without match queryIDs to build all_na_matchDT\n")

all_na_matchDT <- foreach(na_id = noMatch_tads, .combine='rbind') %dopar% {
  id_hicds <- inputGR_DT$HICDS[inputGR_DT$ID == na_id]
  id_exprds <- inputGR_DT$EXPRDS[inputGR_DT$ID == na_id]
  id_tad <- inputGR_DT$TAD[inputGR_DT$ID == na_id]
  stopifnot( paste0(id_hicds, "_", id_exprds, "_", id_tad) == na_id)
  na_matchDT <- data.frame(
    query_id = na_id,
    query_hicds = id_hicds,
    query_exprds = id_exprds,
    queryTAD = id_tad,
    matching_id = NA_character_,
    matching_hicds = NA_character_,
    matching_exprds = NA_character_,
    matchingTAD = NA_character_,
    matchingRatio = NA_real_,
    stringsAsFactors = FALSE
  )
}

all_matchDT <- rbind(all_matchDT_noNA, all_na_matchDT)
stopifnot( length(unique(all_matchDT$query_id)) == nAllTADs)
outFile <- file.path(outFolder, "all_matchDT.Rdata")
save(all_matchDT, file = outFile)
cat(paste0("... written: ", outFile, "\n"))

# outFile <- "INTERSECT_topTADs_ACROSSDS/all_matchDT.Rdata""
load(outFile)

by(all_matchDT, all_matchDT$query_id, function(subDT) {
  subDT
})

# all_matchDT <- all_matchDT[c(1:3,200:205),]

# by(all_matchDT, list(all_matchDT$query_id, all_matchDT$matching_hicds, all_matchDT$matching_exprds), function(subDT) {
#  return(1)
# })
# 
# by(all_matchDT, list(all_matchDT$query_id,all_matchDT$matching_hicds ), function(subDT) {
#   subDT[which.max(subDT$matchingRatio),]
# })

# for a given hicds and exprds -> select the best match TAD
all_bestMatchDT <- do.call(rbind,
        lapply(split(all_matchDT,list(all_matchDT$query_id,all_matchDT$matching_hicds,all_matchDT$matching_exprds ),drop=T), 
          function(subDT) subDT[which.max(subDT$matchingRatio),]))
rownames(all_bestMatchDT) <- NULL

all_bestMatchDT[
  all_bestMatchDT$query_id == "ENCSR312KHQ_SK-MEL-5_40kb_TCGAskcm_lowInf_highInf_chr7_TAD58" &
    all_bestMatchDT$matching_hicds == "K562_40kb" &
    all_bestMatchDT$matching_exprds == "TCGAlaml_wt_mutFLT3" ,
]
all_matchDT[
  all_matchDT$query_id == "ENCSR312KHQ_SK-MEL-5_40kb_TCGAskcm_lowInf_highInf_chr7_TAD58" &
    all_matchDT$matching_hicds == "K562_40kb" &
    all_matchDT$matching_exprds == "TCGAlaml_wt_mutFLT3" ,
  ]

stopifnot(!duplicated(all_bestMatchDT[, c("query_id", "matching_hicds", "matching_exprds")]))

stopifnot( length(unique(all_bestMatchDT$query_id)) == length(unique(na.omit(all_matchDT)$query_id)))

outFile <- file.path(outFolder, "all_bestMatchDT.Rdata")
save(all_bestMatchDT, file = outFile)
cat(paste0("... written: ", outFile, "\n"))

stopifnot( length(unique(all_bestMatchDT$query_id)) == nAllTADs)


# for DEBUG
outFile <- "INTERSECT_topTADs_ACROSSDS/top3/all_matchDT.Rdata"
load(outFile)
all_matchDT[1:5,1:5]

outFile <- "INTERSECT_topTADs_ACROSSDS/top3/all_bestMatchDT.Rdata"
load(outFile)
all_bestMatchDT[1:5,1:5]

# compute the number of datasets in which signif (with or without counting same exrpds)
all_nMatchDT <- aggregate(matching_id ~ query_id, FUN=length, data = all_bestMatchDT)
head(all_nMatchDT)
colnames(all_nMatchDT) <- c("query_id", "all_nMatch")

diffExprds_nMatchDT <- aggregate(matching_id ~ query_id, FUN=length, data = all_bestMatchDT[all_bestMatchDT$query_exprds != all_bestMatchDT$matching_exprds,])
head(diffExprds_nMatchDT)
colnames(diffExprds_nMatchDT) <- c("query_id", "diffExprds_nMatch")

queryID_matchDT <- merge(all_nMatchDT, diffExprds_nMatchDT, by="query_id", all.x = TRUE, all.y=TRUE)
stopifnot(!is.na(queryID_matchDT$all_nMatch))
sum(is.na(queryID_matchDT$diffExprds_nMatch))

stopifnot(!duplicated(queryID_matchDT$query_id))

plot_multiDens(list(
  nMatch_all = queryID_matchDT$all_nMatch,
  nMatch_diffExprds = queryID_matchDT$diffExprds_nMatch),
  plotTit="# datasets with matching signif. TAD", legTxt=NULL, legPos="topright", my_ylab="density", my_xlab="")

plotCex=1.2

plot( 
  x= queryID_matchDT$all_nMatch,
y = queryID_matchDT$diffExprds_nMatch,
ylab = paste0("diffExprds_nMatch"),
xlab = paste0("all_nMatch"),
cex.lab=plotCex,
cex.axis=plotCex,
pch=16,cex=0.7
)

plot_cumMatch <- function(dt, tomatch){
  curr_match <- na.omit(dt[, tomatch])
  xvect <- seq_len(max(curr_match))
  yvect <- sapply(xvect, function(x){
    sum(curr_match >= x)
  })
  plot(x = xvect,
       y = yvect,
       xlab = paste0("# datasets in which matching signif. TAD"), 
       ylab = paste0("# query TAD"),
       type="l")
  mtext(side=3, text=paste0(tomatch))
}
plot_cumMatch(queryID_matchDT, "all_nMatch")
plot_cumMatch(queryID_matchDT, "diffExprds_nMatch")

xvect <- seq_len(max(queryID_matchDT$all_nMatch, na.rm=TRUE))  
yvect <- sapply(xvect, function(x){
  sum(queryID_matchDT$all_nMatch >= x)
})
plot(x = xvect,
     y = yvect,
     xlab = paste0("# datasets in which matching signif. TAD"), 
     ylab = paste0("# query TAD"),
     type="l")
### Problem of dup ??? si un groupe de TADs sont tous best match entre eux -> compté à double dans la courbe ????????



plot_multiDens(list(
  all_matchingRatio = all_matchDT$matchingRatio,
  best_matchingRatio = all_bestMatchDT$matchingRatio),
  plotTit="matchingRatio", legTxt=NULL, legPos="topleft", my_ylab="density", my_xlab="")


source("plot_lolliTAD_funct.R")
top_to_plot = 5

srtd_DT <- queryID_matchDT[order(queryID_matchDT$all_nMatch, decreasing = T),]
head(srtd_DT)
View(srtd_DT)
i=1

# colnames(all_bestMatchDT)
# [1] "query_id"        "query_hicds"     "query_exprds"    "queryTAD"        "matching_id"     "matching_hicds" 
# [7] "matching_exprds" "matchingTAD"     "matchingRatio"


for(i in seq_len(top_to_plot)) {
  
  plot_list <- list()
  
  tad_id <- srtd_DT$query_id[i]
  
  plot_list[[1]] <- plot_lolliTAD_ds(exprds = unique(all_bestMatchDT$query_exprds[all_bestMatchDT$query_id == tad_id]), 
                                  hicds = unique(all_bestMatchDT$query_hicds[all_bestMatchDT$query_id == tad_id]), 
                                  all_TADs = unique(all_bestMatchDT$queryTAD[all_bestMatchDT$query_id == tad_id]))
  
  matching_ids <- all_bestMatchDT$matching_id[all_bestMatchDT$query_id == tad_id]
  stopifnot(!duplicated(matching_ids))
  
  matching_ids <- matching_ids[1:3]
  
  j <- 2
  
  for(matchTAD in matching_ids) {
    plot_list[[j]] <- plot_lolliTAD_ds(exprds = unique(all_bestMatchDT$query_exprds[all_bestMatchDT$query_id == matchTAD]), 
                                       hicds = unique(all_bestMatchDT$query_hicds[all_bestMatchDT$query_id == matchTAD]), 
                                       all_TADs = unique(all_bestMatchDT$queryTAD[all_bestMatchDT$query_id == matchTAD]))
    j <- j+1
  }
  
  all_plots <- do.call(grid.arrange, c(vect_plot,  list(ncol=2, top=textGrob(paste(exprds),
                                                                             gp=gpar(fontsize=20,font=2)))))
  # cat("myHeight =", outHeight, "\n")
  
  outFile <- file.path(outFold, paste0(exprds, "_", hicds, "_", paste0(all_TADs, collapse = "_"), ".", plotType))
  ggsave(filename = outFile, all_plots, width=outWidth, height = outHeight)
  cat("... written: ", outFile, "\n")
  

}



# ######################################################################################
# ######################################################################################

### Pour chaque hicds/exprds -> % de ces signif TADs qui matchent un nombre "x" de signif TADs d'autres datasets

hicds_exprds_DT <- unique(all_bestMatchDT[, c("query_hicds", "query_exprds")])

i=1
for(i in seq_len(nrow(hicds_exprds_DT))) {
  
  signifTADs <- unique(all_bestMatchDT$query_id[all_bestMatchDT$query_hicds == hicds_exprds_DT$query_hicds[i] &
                                           all_bestMatchDT$query_exprds == hicds_exprds_DT$query_exprds[i] 
                                           ])
  
  sapply(signifTADs, function(x) {
    sum(all_bestMatchDT$matching_id == x)
  })
  
}


### Pour chaque hicds/exprds, combien de fois ce dataset est dans le best matching TAD des autres datasets
for(i in seq_len(nrow(hicds_exprds_DT))) {
  
sum(all_bestMatchDT$query_hicds == hicds_exprds_DT$query_hicds[i] &
                                                  all_bestMatchDT$query_exprds == hicds_exprds_DT$query_exprds[i] 
                                                )
}




# ######################################################################################
# ######################################################################################


  
# ######################################################################################
# ######################################################################################
# ######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))



# store the best match in signifTADs_allDS_data, in case needed
# nrow(signifTADs_allDS_data[[file.path(query_hicds, query_exprds)]]$matchDT) == length(signifTADs_allDS_data[[file.path(query_hicds, query_exprds)]]$ids)
# newMatchDT <- 
#   oldMatchDT <- signifTADs_allDS_data[[file.path(query_hicds, query_exprds)]]$matchDT
# oldMatchDT <- oldMatchDT[oldMatchDT$]




# 
# # ensure the matching is done intrachromosomally
# queryIDs_chr <- gsub("(chr.+)_.+", "\\1", queryIDs)
# objectIDs_chr <- gsub("(chr.+)_.+", "\\1", objectIDs)
# stopifnot( queryIDs_chr == objectIDs_chr)
# 
# 
#     query_IR <- IRanges(start = objectTADs_posDT$start,
#                                  width = (objectTADs_posDT$end - objectTADs_posDT$start + 1),
#                                  names=objectTADs_posDT$region)
# 
#     object_allTADs_GR <- GRanges(ranges = query_IR,
#                                    seqnames=gsub("(chr.+)_TAD.+", "\\1", objectTADs_posDT$region))
# 
# 
#     IDoverlap_hits <- findOverlaps(query=query_GR,
#                                subject=object_allTADs_GR)
#     # CHANGED 30.01 _IR to _GR
#     IDoverlaps <- pintersect(query_GR[queryHits(IDoverlap_hits)],
#                            object_allTADs_GR[subjectHits(IDoverlap_hits)])
# 
# 

#   query_IR <- IRanges(start = signifTADs_posDT$start, 
#                            width = (signifTADs_posDT$end - signifTADs_posDT$start + 1), 
#                            names=signifTADs_posDT$region)
#   
#   
#   query_GR <- GRanges(ranges = query_IR,
#                                  seqnames=gsub("(chr.+)_TAD.+", "\\1", signifTADs_posDT$region))
#   

# ######################################################################################
# cat(paste0("... written: ", outFile, "\n"))
# cat(paste0("... written: ", listFile, "\n"))



# hicds="ENCSR444WCZ_A549_40kb"
# all_topTADsGenes <- foreach(hicds = all_hicds) %dopar% {
#   
#   if(hicds == "pipelineConsensus") {
#     dsPipOutDir <- file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom", "OUTPUT_FOLDER",  exprds)
#     
#     # file with assignment from entrez to all regions
#     gene2tadDT_file <- paste0(setDir, 
#                               "/mnt/ed4/marie/gene_data_final/consensus_TopDom_covThresh_r0.6_t80000_v0_w-1_final/genes2tad/all_genes_positions.txt")
#     TADposDT_file <- paste0(setDir, 
#                               "/mnt/ed4/marie/gene_data_final/consensus_TopDom_covThresh_r0.6_t80000_v0_w-1_final/genes2tad/all_assigned_regions.txt") 
#     
#   } else {
#     stopifnot(dir.exists(hicds))
#     dsPipOutDir <- file.path("PIPELINE", "OUTPUT_FOLDER", hicds, exprds)
#     gene2tadDT_file <- file.path(hicds, "genes2tad/all_genes_positions.txt")
#     TADposDT_file <- file.path(hicds, "genes2tad/all_assigned_regions.txt")
#   }
#   stopifnot(file.exists(gene2tadDT_file))
#   g2tDT <- read.delim(gene2tadDT_file, stringsAsFactors = FALSE, header=F, col.names=c("entrezID", "chromo", "start", "end", "region"))
#   
#   stopifnot(file.exists(TADposDT_file))
#   TADposDT <- read.delim(TADposDT_file, stringsAsFactors = FALSE, header=F, col.names=c("chromo","region", "start", "end"))
#   stopifnot(is.numeric(TADposDT$start))
#   TADposDT <- TADposDT[grepl("_TAD", TADposDT$region),]
#   stopifnot(nrow(TADposDT) > 0)
#   
#   stopifnot(g2tDT$entrezID %in% entrezDT$entrezID)
#   g2t2sDT <- merge(g2tDT[, c("entrezID", "region")], entrezDT[,c("entrezID", "symbol")], by="entrezID", all.x = TRUE, all.y = FALSE )
#   head(g2t2sDT)
#   stopifnot(!is.na(g2t2sDT))
#   g2t2sDT$entrezID <- as.character(g2t2sDT$entrezID)
#   g2t2sDT$region <- as.character(g2t2sDT$region)
#   g2t2sDT$symbol <- as.character(g2t2sDT$symbol)
#   
#   # subset the genes used in the pipeline
#   stopifnot(dir.exists(dsPipOutDir))
#   stopifnot(dir.exists(file.path(dsPipOutDir, script0_name)))
#   # geneList: values are entrez, names can be ENSEMBL, entrez, etc.
#   geneList_file <- file.path(dsPipOutDir, script0_name, "pipeline_geneList.Rdata")
#   stopifnot(file.exists(geneList_file))
#   geneList <- eval(parse(text = load(geneList_file)))
#   stopifnot(geneList %in% g2tDT$entrezID)
#   g2tDT <- g2tDT[g2tDT$entrezID %in% geneList,]
#   
#   stopifnot(dir.exists(file.path(dsPipOutDir, script11_name)))
#   tad_pvalFile <- file.path(dsPipOutDir, script11_name, "emp_pval_combined.Rdata")
#   stopifnot(file.exists(tad_pvalFile))
#   tad_pvals <- eval(parse(text = load(tad_pvalFile)))
#   
#   adj_tad_pvals <- sort(p.adjust(tad_pvals, method="BH"))
#   
#   if(topThresh >= 1) {
#     topThresh <- min(c(topThresh, length(adj_tad_pvals)))
#     pvalThresh <- as.numeric(adj_tad_pvals[topThresh])
#     stopifnot(!is.na(pvalThresh))
#     stopifnot(pvalThresh <= 1)
#   } else {
#     pvalThresh <- topThresh
#   }
#   top_pvals <- adj_tad_pvals[adj_tad_pvals <= pvalThresh]
#   stopifnot(!is.na(top_pvals))
#   
#   if(topThresh >= 1) stopifnot(length(unique(top_pvals)) <= topThresh)
#   
#   top_tads <- names(top_pvals)
#   
#   tad_ranks <- rank(top_pvals, ties="min")
#   stopifnot(top_tads %in% names(tad_ranks))
#   tad_ranks <- tad_ranks[top_tads]
#   
#   topTADsDT <- data.frame(
#     top_tads = top_tads,
#     tad_rank = tad_ranks,
#     top_pvals = as.numeric(top_pvals),
#     stringsAsFactors = FALSE
#   )
#   head(topTADsDT)
#   
#   stopifnot(top_tads %in% g2t2sDT$region)
#   topGenesDT <- g2t2sDT[g2t2sDT$region %in% top_tads,]
#   head(topGenesDT)
#   
#   colnames(topGenesDT)[colnames(topGenesDT) == "region"] <- "top_tads"
#   
#   topTADs_genesDT <- merge(topTADsDT, topGenesDT[, c("top_tads", "symbol")], 
#                            by = "top_tads", all.x=TRUE, all.y=TRUE)
#   
#   topTADs_genesDT$top_tads <- as.character(topTADs_genesDT$top_tads)
#   topTADs_genesDT$symbol <- as.character(topTADs_genesDT$symbol)
#   topTADs_genesDT <- topTADs_genesDT[order(topTADs_genesDT$top_pvals, topTADs_genesDT$top_tads, topTADs_genesDT$symbol),]
#   head(topTADs_genesDT)
#   
#   list(topTADs_genesDT=topTADs_genesDT,
#        TADposDT=TADposDT,
#        gene2tadDT = g2tDT)
# }
# names(all_topTADsGenes) <- all_hicds
# str(all_topTADsGenes)
# 
# ########################################################################### LOAD SAVED DATA
# # outFolder <- "INTERSECT_topTADs"
# # outFile <- file.path(outFolder, "all_topTADsGenes.Rdata")
# # save(all_topTADsGenes, file = outFile)
#       # outFile = "INTERSECT_topTADs/all_topTADsGenes.Rdata"
#       # load(outFile)
#       # load("INTERSECT_topTADs/all_topTADsGenes.Rdata")
#       # str(all_topTADsGenes)
# 
# # names(all_topTADsGenes[[1]])
# # [1] "topTADs_genesDT" "TADposDT"        "gene2tadDT" 
# 
# i=1
# j=2
# 
# all_maxIntersectOverlap_matchDT <- foreach(i = seq_along(all_topTADsGenes), .combine="rbind") %dopar% {
#   
#   ref_ds <- names(all_topTADsGenes)[i]
#   
#   cat("... start with ref_ds\t=\t", ref_ds, "\n")
#   
#   signifTADs <- unique(as.character(all_topTADsGenes[[i]][["topTADs_genesDT"]][, "top_tads"]))
#   
#   signifTADs_posDT <- all_topTADsGenes[[i]][["TADposDT"]]
#   signifTADs_posDT$region <- as.character(signifTADs_posDT$region)
#   stopifnot(signifTADs %in% signifTADs_posDT$region)
#   signifTADs_posDT <- signifTADs_posDT[signifTADs_posDT$region %in% signifTADs,]
#   
#   ref_g2tDT <- all_topTADsGenes[[i]][["gene2tadDT"]]
#   
#   query_IR <- IRanges(start = signifTADs_posDT$start, 
#                            width = (signifTADs_posDT$end - signifTADs_posDT$start + 1), 
#                            names=signifTADs_posDT$region)
#   
#   
#   query_GR <- GRanges(ranges = query_IR,
#                                  seqnames=gsub("(chr.+)_TAD.+", "\\1", signifTADs_posDT$region))
#   
#   
#   ref_maxIntersectOverlap_matchDT <- foreach(j = seq_along(all_topTADsGenes)[-i], .combine="rbind") %dopar% {
#     
#     obj_ds <- names(all_topTADsGenes)[j]
#     
#     cat("... start with obj_ds\t=\t", obj_ds, "\n")
#     
#     objectTADs_posDT <- all_topTADsGenes[[j]][["TADposDT"]]
#     objectTADs_posDT$region <- as.character(objectTADs_posDT$region)
#     
#     obj_g2tDT <- all_topTADsGenes[[j]][["gene2tadDT"]]
#     obj_g2tDT$region <- as.character(obj_g2tDT$region)
#     # added 30.01: filter the TAD to the ones in g2t -> I want only those with genes inside
#     objectTADs_posDT <- objectTADs_posDT[objectTADs_posDT$region %in% obj_g2tDT$region, ]
#     stopifnot(nrow(objectTADs_posDT) > 0)
#     
#     query_IR <- IRanges(start = objectTADs_posDT$start, 
#                                  width = (objectTADs_posDT$end - objectTADs_posDT$start + 1), 
#                                  names=objectTADs_posDT$region)
#     
#     object_allTADs_GR <- GRanges(ranges = query_IR,
#                                    seqnames=gsub("(chr.+)_TAD.+", "\\1", objectTADs_posDT$region))
#     
#     
#     IDoverlap_hits <- findOverlaps(query=query_GR, 
#                                subject=object_allTADs_GR)
#     # CHANGED 30.01 _IR to _GR
#     IDoverlaps <- pintersect(query_GR[queryHits(IDoverlap_hits)], 
#                            object_allTADs_GR[subjectHits(IDoverlap_hits)])
#     
#     queryIDs <- names(query_IR[queryHits(IDoverlap_hits)])
#     objectIDs <- names(query_IR[subjectHits(IDoverlap_hits)])
#     stopifnot( length(queryIDs) == length(objectIDs) )
#     
#     # ensure the matching is done intrachromosomally
#     queryIDs_chr <- gsub("(chr.+)_.+", "\\1", queryIDs)
#     objectIDs_chr <- gsub("(chr.+)_.+", "\\1", objectIDs)
#     stopifnot( queryIDs_chr == objectIDs_chr)
#     
#     # retrieve the number of common genes for each match
#     # take the reference TADs that have a match
#     queryTADs_withMatch <- unique(queryIDs)
#     
#     txt <- paste0("... signif. queryIDs with overlapping objectIDs:\t")
#     printAndLog(txt, logFile)
#     txt <- paste0(length(queryTADs_withMatch), "/", length(query_GR), "\n")
#     printAndLog(txt, logFile)
#     
#     cat("... retrieve genes for the TADs that have matching object TADs\n")
#     
#     
#     # for each TAD, 1) retrieve the genes that belong to it
#     queryTADs_withMatch_genes <- lapply(queryTADs_withMatch, function(curr_tad){
#       stopifnot( curr_tad %in% ref_g2tDT$region )
#       as.character(ref_g2tDT$entrezID[ref_g2tDT$region == curr_tad])
#     })
#     names(queryTADs_withMatch_genes) <- queryTADs_withMatch
#     
#     # $chr9_TAD191
#     # [1] "9858"   "138162" "51116"  "3933"   "29991"  "5047"   "402381" "57582"  "157922" "10422"  "138151" "90120"  "169714"
#     
#     # queryTADs_withMatch <- queryTADs_withMatch[1:3]
#     
#     # for each TAD, 2) retrieve the genes of that TAD that matches with it
#     
#     cat("... retrieve the genes of the matching TADs\n")
#     
#     queryTADs_withMatch_matchingObjectTADs <- lapply(queryTADs_withMatch, function(curr_tad){
#           
#           ref_genes <- queryTADs_withMatch_genes[[curr_tad]]
#           
#           matching_objectTADs <- objectIDs[queryIDs == curr_tad]
#           stopifnot(length(matching_objectTADs) > 0)
#           
#           matching_objectTADs_genes <- lapply(matching_objectTADs, function(obj_tad){
#             stopifnot( obj_tad %in% obj_g2tDT$region )
#             genes_in_matchingTADs <- as.character(obj_g2tDT$entrezID[obj_g2tDT$region == obj_tad])
#             nbr_intersectGenes_in_matchingTADs <- sum(genes_in_matchingTADs %in% ref_genes)
#             list(matchingTADs_genes = genes_in_matchingTADs,
#                  matchingTADs_nIntersect = nbr_intersectGenes_in_matchingTADs)
#           })
#           names(matching_objectTADs_genes) <- matching_objectTADs
#           matching_objectTADs_genes
#         
#           # list(matching_objectTADs=matching_objectTADs,
#           #      matching_objectTADs_genes=matching_objectTADs_genes)
#     })
#     names(queryTADs_withMatch_matchingObjectTADs) <- queryTADs_withMatch
#     
#     # queryTADs_withMatch_matchingObjectTADs[["chr9_TAD191"]]
#     # $chr9_TAD178
#     # $chr9_TAD178$matchingTADs_genes
#     # [1] "1289"  "2220"  "2219"  "10439"
#     # 
#     # $chr9_TAD178$matchingTADs_nIntersect
#     # [1] 0
#     
#     ### VERSION 1 -> MATCHING TAD IS THE ONE WITH MOST INTERSECT GENES
#     
#     # e.g. for the query chr10_TAD1 -> match 2 objectIDs
#     # > queryTADs_withMatch_matchingObjectTADs[["chr10_TAD1"]]
#     # $chr10_TAD1
#     # $chr10_TAD1$matchingTADs_genes
#     # [1] "347688"    "439945"    "10771"     "100421369"
#     # 
#     # $chr10_TAD1$matchingTADs_nIntersect
#     # [1] 4
#     # 
#     # 
#     # $chr10_TAD2
#     # $chr10_TAD2$matchingTADs_genes
#     # [1] "22982"     "100847086" "414235"   
#     # 
#     # $chr10_TAD2$matchingTADs_nIntersect
#     # [1] 3
#     
#     # > queryTADs_withMatch_matchingObjectTADs[["chr11_TAD130"]]
#     # $chr11_TAD108
#     # $chr11_TAD108$matchingTADs_genes
#     # [1] "246330"    "101928069" "10072"     "582"       "254359"    "89"       
#     # [7] "8722"      "55231"     "9973"      "10432"     "100526737" "5936"     
#     # [13] "83759"    
#     # 
#     # $chr11_TAD108$matchingTADs_nIntersect
#     # [1] 13
#     # 
#     # 
#     # $chr11_TAD109
#     # $chr11_TAD109$matchingTADs_genes
#     # [1] "6712"      "79703"     "100422299" "100462788"
#     # 
#     # $chr11_TAD109$matchingTADs_nIntersect
#     # [1] 4
#     
#     # Filter to retain only the highest intersect
#     
#     cat("... filter retain highest gene intersect\n")
#     
#     queryTADs_withMatch_maxInterObjectTADs <- lapply(queryTADs_withMatch_matchingObjectTADs, function(curr_tad){
#       tmp_nIntersect <- unlist(lapply(curr_tad, function(x) as.numeric(x[["matchingTADs_nIntersect"]])))
#       stopifnot( ! is.na(tmp_nIntersect) )
#       obj_maxIntersect_names <- which.max(tmp_nIntersect)
#       obj_maxIntersect_names
#     })
#     stopifnot(names(queryTADs_withMatch_maxInterObjectTADs) == queryTADs_withMatch)
#     
#     maxIntersectMatchingTAD <- as.character(unlist(
#       lapply(queryTADs_withMatch_maxInterObjectTADs, function(x) names(x))))
#     
#     nIntersectMatchingTAD <- unlist(lapply(queryTADs_withMatch_matchingObjectTADs, function(curr_tad){
#       tmp_nIntersect <- unlist(lapply(curr_tad, function(x) as.numeric(x[["matchingTADs_nIntersect"]])))
#       stopifnot( ! is.na(tmp_nIntersect) )
#       obj_maxIntersect_nbr <- max(tmp_nIntersect)
#       obj_maxIntersect_nbr
#     }))
# 
#     stopifnot(!is.na(nIntersectMatchingTAD))
#     stopifnot(is.numeric(nIntersectMatchingTAD))
#     
#     stopifnot( length(maxIntersectMatchingTAD) == length(queryTADs_withMatch) )
#     
#     query_nbrGenes <- sapply(queryTADs_withMatch, function(x) 
#                           sum(ref_g2tDT$region == x))
#     
#     object_nbrGenes <- sapply(maxIntersectMatchingTAD, function(x) 
#       sum(obj_g2tDT$region == x))
#     
#     query_matching_ratio <- nIntersectMatchingTAD/query_nbrGenes
#     stopifnot(query_matching_ratio >= 0 & query_matching_ratio <= 1 )
#     
#     maxIntersect_matchDT <- data.frame(
#       
#       expr_ds = exprds,
#       
#       query_ds = ref_ds,
#       matching_ds = obj_ds,
#       
#       queryTAD = queryTADs_withMatch,
#       query_nbrGenes = query_nbrGenes,
#       
#       matchingTAD_maxIntersect = maxIntersectMatchingTAD,
#       matching_nbrGenes = object_nbrGenes,
#         
#       query_matching_intersectGenes = nIntersectMatchingTAD,
#       query_matching_ratio = query_matching_ratio,
#       
#       
#       stringsAsFactors = FALSE
#     )
#     stopifnot(maxIntersect_matchDT$query_matching_intersectGenes <= maxIntersect_matchDT$matching_nbrGenes)
#     stopifnot(maxIntersect_matchDT$query_matching_intersectGenes <= maxIntersect_matchDT$query_nbrGenes)
#     
#     maxIntersect_matchDT$queryTAD <- as.character(maxIntersect_matchDT$queryTAD)
#     maxIntersect_matchDT <- maxIntersect_matchDT[order(maxIntersect_matchDT$queryTAD),]
#     
#                 # outFile <- file.path(outFolder, paste0(exprds, "_", ref_ds, "_", obj_ds, "_maxIntersect_matchDT.txt"))
#                 # cat("... writte table with max gene intersect\n")
#                 # write.table(maxIntersect_matchDT, col.names=TRUE, row.names=FALSE, sep="\t", quote=F, file = outFile)
#                 # cat(paste0("... written: ", outFile, "\n"))
#         
#     ### VERSION 2 -> MATCHING TAD IS THE ONE WITH HIGHEST OVERLAP
#         
#     # IDoverlap_hits <- findOverlaps(query=query_GR, 
#     #                                 subject=object_allTADs_GR)
#     # 
#     # IDoverlaps <- pintersect(query_GR[queryHits(IDoverlap_hits)], 
#     #                           object_allTADs_GR[subjectHits(IDoverlap_hits)])
#     percentTADoverlap <- width(IDoverlaps)/width(query_GR[queryHits(IDoverlap_hits)])
#     stopifnot( percentTADoverlap > 0 & percentTADoverlap <= 1 )
#     
# 
#     overlap_matchDT <- data.frame(
#       
#       expr_ds = exprds,
#       
#       query_ds = ref_ds,
#       matching_ds = obj_ds,
#       
#       queryTAD = names(query_GR[queryHits(IDoverlap_hits)]),
#       query_size = width(query_GR[queryHits(IDoverlap_hits)]),
#       
#       matchingTAD = names(object_allTADs_GR[subjectHits(IDoverlap_hits)]),
#       matching_size = width(object_allTADs_GR[subjectHits(IDoverlap_hits)]),
#       
#       query_matching_overlap = percentTADoverlap,
#       stringsAsFactors = FALSE
#     )
#     
# 
# 
#         
#     maxOverlap_matchDT <- do.call(rbind, by(data = overlap_matchDT, INDICES = overlap_matchDT$queryTAD, FUN=function(x) {
#       x[which.max(x$query_matching_overlap),]
#     }))
# 
#     colnames(maxOverlap_matchDT)[  colnames(maxOverlap_matchDT) == "matchingTAD"] <- "matchingTAD_maxOverlap"
#     
#     maxOverlap_matchDT$queryTAD <- as.character(maxOverlap_matchDT$queryTAD)
#     maxOverlap_matchDT <- maxOverlap_matchDT[order(maxOverlap_matchDT$queryTAD),]
#     
#                 # outFile <- file.path(outFolder, paste0(exprds, "_", ref_ds, "_", obj_ds, "_maxOverlap_matchDT.txt"))
#                 # cat("... writte table with max overlap\n")
#                 # write.table(maxOverlap_matchDT, col.names=TRUE, row.names=FALSE, sep="\t", quote=F, file = outFile)
#                 # cat(paste0("... written: ", outFile, "\n"))
#     
#     
#     # list(
#     #   maxIntersect = maxIntersect_matchDT,
#     #   maxOverlap = maxOverlap_matchDT
#     # )
#     #break
#     
#     # expr_ds	query_ds	matching_ds	queryTAD	query_nbrGenes	matchingTAD_maxIntersect	matching_nbrGenes	query_matching_intersectGenes
#     # expr_ds	query_ds	matching_ds	queryTAD	query_size	matchingTAD_maxOverlap	matching_size	query_matching_overlap
#     
#     maxIntersectOverlap_matchDT <- merge(maxIntersect_matchDT, maxOverlap_matchDT, by=c("expr_ds", "query_ds", "matching_ds", "queryTAD"))
#     maxIntersectOverlap_matchDT
#     
#   } # end iterating over the j [object dataset] -> ref_maxIntersectOverlap_matchDT
#   #break
#   ref_maxIntersectOverlap_matchDT
# } # end iterating over the i [reference dataset] -> all_maxIntersectOverlap_matchDT
# 
# outFile <- file.path(outFolder, paste0(exprds, "_", ref_ds, "_", obj_ds, "_top", topThresh, "_all_maxIntersectOverlap_matchDT.Rdata"))
# cat("... writte table with max overlap and intersect\n")
# save(all_maxIntersectOverlap_matchDT,file = outFile)
# cat(paste0("... written: ", outFile, "\n"))
# 
# 
# outDT <- all_maxIntersectOverlap_matchDT
# 
# # outDT <- outDT[order(outDT[, "query_matching_ratio"]),]
# 
# 
# outDT[, "query_matching_ratio"] <- round(outDT[, "query_matching_ratio"], 4)
# outDT[, "query_matching_overlap"] <- round(outDT[,"query_matching_overlap"], 4)
# 
# 
# outFile <- file.path(outFolder, paste0(exprds, "_", ref_ds, "_", obj_ds, "_top", topThresh, "_all_maxIntersectOverlap_matchDT.txt"))
# cat("... writte table with max overlap and intersect\n")
# write.table(outDT, col.names=TRUE, row.names=FALSE, sep="\t", quote=F, file = outFile)
# cat(paste0("... written: ", outFile, "\n"))
# 
# shortDT <- outDT 
# shortDT <- shortDT[, c("expr_ds", "query_ds", "matching_ds", "queryTAD", "matchingTAD_maxIntersect", "query_matching_ratio")]
# shortDT <- shortDT[order(shortDT[, "query_matching_ratio"], decreasing = TRUE),]
# outFile <- file.path(outFolder, paste0(exprds, "_", ref_ds, "_", obj_ds, "_top", topThresh, "_all_maxIntersectOverlap_matchDT_shortDT.txt"))
# cat("... writte table with max overlap and intersect\n")
# write.table(shortDT, col.names=TRUE, row.names=FALSE, sep="\t", quote=F, file = outFile)
# cat(paste0("... written: ", outFile, "\n"))







