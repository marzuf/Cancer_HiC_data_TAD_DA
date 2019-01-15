
# Rscript check_step17_files.R MCF-7_40kb TCGAbrca_lum_bas

hicdt="MCF-7_40kb"
exprdt="TCGAbrca_lum_bas"

setDir <- "~/media/electron"
setDir <- ""
setwd(file.path(setDir, "/mnt/etemp/marie/Cancer_HiC_data_TAD_DA"))

args <- commandArgs(trailingOnly = TRUE)
hicdt <- args[1]
exprdt <- args[2]

all_ratiosFile <- file.path("PIPELINE", "OUTPUT_FOLDER", 
                            hicdt, exprdt, 
                            "170revision2EZH2_score_auc_pval_permGenes/auc_ratios.Rdata")
stopifnot(file.exists(all_ratiosFile))

all_ratios <- eval(parse(text=load(all_ratiosFile)))

all_ratios
stopifnot(names(all_ratios) == c("ratioDown_auc_permGenes","ratioDown_auc_permGenes","prodSignedRatio_auc_permGenes","prodSignedRatio_auc_permGenes"))

cat("*** DONE - ok\n")