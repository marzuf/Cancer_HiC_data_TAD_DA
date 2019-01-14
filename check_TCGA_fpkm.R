
# Rscript check_TCGA_fpkm.R MCF-7_40kb TCGAbrca_lum_bas

setDir <- "~/media/electron"
setDir <- ""
setwd(file.path(setDir, "/mnt/etemp/marie/Cancer_HiC_data_TAD_DA"))

args <- commandArgs(trailingOnly = TRUE)
hicdt <- args[1]
exprdt <- args[2]

# in setting file:
# UPDATE 07.12.2018: for RSEM data, the "analog" FPKM file is provided separately (built in prepData)
rna_fpkmDT_file <- file.path(setDir, "/mnt/ed4/marie/other_datasets", exprdt, "fpkmDT.Rdata")
settingDT <- eval(parse(text = load(rna_fpkmDT_file)))
dim(settingDT)

pip_fpkm_file <- file.path("PIPELINE", "OUTPUT_FOLDER", hicdt, exprdt, "0_prepGeneData", "rna_fpkmDT.Rdata")
pipDT <- eval(parse(text = load(pip_fpkm_file)))
dim(pipDT)

commonCols <- intersect(colnames(pipDT), colnames(settingDT))

settingDT <- settingDT[, commonCols]
pipDT <-  pipDT[, commonCols]


commonRows <- intersect(rownames(pipDT), rownames(settingDT))
settingDT <- settingDT[commonRows,]
pipDT <-  pipDT[commonRows,]

all.equal(settingDT, pipDT)

stopifnot(all.equal(settingDT, pipDT))

cat("*** DONE - ok\n")