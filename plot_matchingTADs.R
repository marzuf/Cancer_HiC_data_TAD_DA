# Rscript intersect_topTADs_acrossDS.R 3

startTime <- Sys.time()

cat("> START plot_matchingTADs.R \n")

# Rscript plot_matchingTADs.R NCI-H460_40kb_TCGAlusc_norm_lusc_chr1_TAD232 INTERSECT_topTADs_ACROSSDS/top3/all_matchDT.Rdata
# Rscript plot_matchingTADs.R NCI-H460_40kb_TCGAlusc_norm_lusc_chr1_TAD232 3

SSHFS <- FALSE
setDir <- ifelse(SSHFS, "~/media/electron", "")

plotType <- "svg"

require(ggplot2)
source( file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg"), "set_dataset_colors.R"))
head(score_DT)
source( file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_11_18"), "analysis_utils.R"))
source( file.path("colors_utils.R"))
dataset_proc_colors <- setNames(score_DT$proc_col, score_DT$dataset)
length(dataset_proc_colors)

source("utils_fct.R")
source("plot_lolliTAD_funct.R")

outFolder <- file.path("PLOT_MATCHINGTADs")
dir.create(outFolder)

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) > 0)
tad_id <- args[1]
dt_arg <- as.numeric(args[2])
if(is.na(dt_arg)) {
  matchDT_file <- args[2]
} else {
  matchDT_file <- file.path("INTERSECT_topTADs_ACROSSDS", paste0("top", dt_arg), "all_matchDT.Rdata")
}
stopifnot(file.exists(matchDT_file))
all_matchDT <- eval(parse(text = load(matchDT_file)))

all_bestMatchDT <- do.call(rbind,
                           lapply(split(all_matchDT,list(all_matchDT$query_id,all_matchDT$matching_hicds,all_matchDT$matching_exprds ),drop=T), 
                                  function(subDT) subDT[which.max(subDT$matchingRatio),]))
rownames(all_bestMatchDT) <- NULL

stopifnot(tad_id %in% all_bestMatchDT$query_id)

plot_list <- list()
  


plot_list[[1]] <- plot_lolliTAD_ds(exprds = unique(all_bestMatchDT$query_exprds[all_bestMatchDT$query_id == tad_id]), 
                                 hicds = unique(all_bestMatchDT$query_hicds[all_bestMatchDT$query_id == tad_id]), 
                                 all_TADs = unique(all_bestMatchDT$queryTAD[all_bestMatchDT$query_id == tad_id]))

matching_ids <- all_bestMatchDT$matching_id[all_bestMatchDT$query_id == tad_id]
stopifnot(!duplicated(matching_ids))


j <- 2
for(matchTAD in matching_ids) {
  plot_list[[j]] <- plot_lolliTAD_ds(exprds = unique(all_bestMatchDT$query_exprds[all_bestMatchDT$query_id == matchTAD]), 
                                     hicds = unique(all_bestMatchDT$query_hicds[all_bestMatchDT$query_id == matchTAD]), 
                                     all_TADs = unique(all_bestMatchDT$queryTAD[all_bestMatchDT$query_id == matchTAD]))
  j <- j+1
}

all_plots <- do.call(grid.arrange, c(plot_list,  list(ncol=2, top=textGrob(paste(tad_id),
                                                                         gp=gpar(fontsize=20,font=2)))))
# cat("myHeight =", outHeight, "\n")
outWidth <- 20
outHeight <- min(c(7 * length(plot_list)/2, 49))

outFile <- file.path(outFolder, paste0(tad_id, "_acrossDsMatchingTADs_lolliplots.", plotType))
ggsave(filename = outFile, all_plots, width=outWidth, height = outHeight)
cat("... written: ", outFile, "\n")



# ######################################################################################
# ######################################################################################
# ######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))
