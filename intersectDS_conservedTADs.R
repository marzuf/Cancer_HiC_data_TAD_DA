# Rscript intersectDS_conservedTADs.R

startTime <- Sys.time()

source("utils_fct.R")

plotCex <- 1.2

plotType <- "png"
myHeight <- ifelse(plotType == "png", 400, 7)
myWidth <- ifelse(plotType == "png", 400, 7)


outFold <- file.path("INTERSECTDS_CONSERVEDTADS")
dir.create(outFold, recursive = TRUE)

inFold <- file.path("INTERSECT_topTADs_ACROSSDS", "top1")
stopifnot(dir.exists(inFold))

cat("load signifTAD\n")
load(file.path(inFold, "signifTADs_allDS_data.Rdata"))
nTADs <- unlist(lapply(signifTADs_allDS_data, function(x) {
  nrow(x[["posDT"]])
}))
totTADs <- sum(nTADs)

head(signifTADs_allDS_data[[1]][["ids"]])
head(signifTADs_allDS_data[[1]][["idDT"]])
head(signifTADs_allDS_data[[1]][["geneDT"]])
head(signifTADs_allDS_data[[1]][["signifDT"]])
head(signifTADs_allDS_data[[1]][["posDT"]])
head(signifTADs_allDS_data[[1]][["matchDT"]])

cat("load all_match\n")
load(file.path(inFold, "all_matchDT.Rdata"))
stopifnot( totTADs == length(unique(all_matchDT$query_id)))
head(all_matchDT)

cat("load bestMatch\n")
load(file.path(inFold, "all_bestMatchDT.Rdata"))
# stopifnot( totTADs == length(unique(all_bestMatchDT$query_id))) # not TRUE, NA removed
head(all_bestMatchDT)

# compute the number of datasets in which signif (with or without counting same exrpds)
all_nMatchDT <- aggregate(matching_id ~ query_id, FUN=length, 
                          data = all_bestMatchDT)
head(all_nMatchDT)
colnames(all_nMatchDT) <- c("query_id", "all_nMatch")

outFile <- file.path(outFold, "all_nMatchDT.Rdata" )
dir.create(dirname(outFile), recursive=TRUE)
save(all_nMatchDT, file=outFile)

all_nMatchDT <- all_nMatchDT[order(all_nMatchDT$all_nMatch, 
                                   decreasing = TRUE),]
head(all_nMatchDT)
tail(all_nMatchDT)
plot(density(all_nMatchDT$all_nMatch))

nGenesByTADs <- unlist(unname(lapply(signifTADs_allDS_data, function(x) {
  tmpDT <- x[["geneDT"]]
  setNames(as.numeric(table(tmpDT$ID)), 
           as.character(names(table(tmpDT$ID))))
})))
nGenesByTADs_DT <- data.frame(
  query_id = names(nGenesByTADs),
  nGenes = nGenesByTADs,
  stringsAsFactors = FALSE
)
rownames(nGenesByTADs_DT) <- NULL

sizeTADs <- unlist(unname(lapply(signifTADs_allDS_data, function(x) {
  tmpDT <- x[["posDT"]]
  setNames(tmpDT$END-tmpDT$START, tmpDT$ID)
})))
sizeTADs_DT <- data.frame(
  query_id = names(sizeTADs),
  TADsize = sizeTADs,
  stringsAsFactors = FALSE
)
rownames(sizeTADs_DT) <- NULL

TADpvals <- unlist(unname(lapply(signifTADs_allDS_data, function(x) {
  tmpDT <- x[["signifDT"]]
  setNames(tmpDT$PVAL, tmpDT$ID)
})))
TADpvals_DT <- data.frame(
  query_id = names(TADpvals),
  TADpval = TADpvals,
  stringsAsFactors = FALSE
)
rownames(TADpvals_DT) <- NULL

features_DT <- merge(merge(merge(all_nMatchDT, nGenesByTADs_DT, by="query_id"),
                     sizeTADs_DT, by = "query_id"), TADpvals_DT, by = "query_id")

stopifnot(nrow(features_DT) > 0)

var_to_plot="nGenes"
for(var_to_plot in c("nGenes", "TADsize", "TADpval")) {
  
  xvar <- var_to_plot
  yvar <- "all_nMatch"
  
  myxlab <- paste0(var_to_plot)
  myylab <- paste0("# datasets matching")
  
  outFile <- file.path(outFold, paste0(yvar, "_vs_", xvar, ".", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  densplot(
    x = features_DT[, xvar],
    y = features_DT[, yvar],
    xlab = myxlab,
    ylab = myylab,
    cex.lab = plotCex,
    cex.axis = plotCex,
    cex = 0.7,
    pch = 16
  )
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  outFile <- file.path(outFold, paste0(yvar, "_vs_", xvar, "_log10.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  densplot(
    x = log10(features_DT[, xvar]),
    y = features_DT[, yvar],
    xlab = paste0(myxlab, " [log10]"),
    ylab = myylab,
    cex.lab = plotCex,
    cex.axis = plotCex,
    cex = 0.7,
    pch = 16
  )
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
}


cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))
