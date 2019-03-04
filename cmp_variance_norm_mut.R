
SSHFS <- TRUE
setDir <- ifelse(SSHFS, "~/media/electron", "")

varianceDT_file <- "GENE_VARIANCE/LOG2FPKM/all_ds_geneVarDT.Rdata"
load(varianceDT_file)

dataset <- "luad"
dataset <- "lihc"

ds_pipFolder <- file.path("PIPELINE", "OUTPUT_FOLDER")
ds_settingFolder <- file.path("PIPELINE", "INPUT_FILES")

if(dataset == "luad"){
  exprds1 <- "TCGAluad_norm_luad"
  exprds2 <- "TCGAluad_mutKRAS_mutEGFR"
  exprds3 <- "TCGAluad_nonsmoker_smoker"
  exprds4 <- "TCGAluad_wt_mutKRAS"
  # otherDS <- c(exprds2, exprds3, exprds4)
  otherDS <- c(exprds2, exprds4)
  hicds <- "NCI-H460_40kb"
  
}else if(dataset == "lihc"){
  exprds1 <- "TCGAlihc_norm_lihc"  
  exprds2 <- "TCGAlihc_wt_mutCTNNB1"
  otherDS <- c(exprds2)
  hicds <- "GSE105381_HepG2_40kb"
}
normDS <- exprds1

######################## VARIANCE COMPUTED

varName <- "meanMostVar"
cat(paste0("For: ", varName, "\n"))
normDS_var <- all_ds_geneVarDT[all_ds_geneVarDT$hicds == hicds & all_ds_geneVarDT$exprds == normDS, varName]
cat(paste0("... ", normDS, ":\t", round(normDS_var, 2 ), "\n"))
for(exprds in otherDS) {
  epxrDS_var <- all_ds_geneVarDT[all_ds_geneVarDT$hicds == hicds & all_ds_geneVarDT$exprds == exprds, varName]
  cat(paste0("... ", exprds, ":\t", round(epxrDS_var, 2), "\n"))  
}


######################## CHECK THE SAMPLES

fpkm_file <- file.path(ds_pipFolder, hicds, normDS, "0_prepGeneData", paste0("rna_", "fpkm", "DT.Rdata"))
stopifnot(file.exists(fpkm_file))

setting_file <- file.path(ds_settingFolder, hicds, paste0("run_settings_", normDS, ".R"))
stopifnot(file.exists(setting_file))
source(setting_file)

norm_samp1 <- eval(parse(text=load(file.path(setDir, sample1_file))))
norm_samp2 <- eval(parse(text=load(file.path(setDir, sample2_file))))

stopifnot(cond1 == "norm" )  
norm_cond2 <- cond2

exprds <- otherDS[1]
for(exprds in otherDS) {
  # fpkm_file <- file.path(ds_pipFolder, hicds, exprds, "0_prepGeneData", paste0("rna_", "fpkm", "DT.Rdata"))
  # stopifnot(file.exists(fpkm_file))
  
  setting_file <- file.path(ds_settingFolder, hicds, paste0("run_settings_", exprds, ".R"))
  stopifnot(file.exists(setting_file))
  source(setting_file)
  
  samp1 <- eval(parse(text=load(file.path(setDir, sample1_file))))
  samp2 <- eval(parse(text=load(file.path(setDir, sample2_file))))
  
  stopifnot(samp1 %in% norm_samp2)
  cat("-> all ", cond1, " in ", norm_cond2, ":\tok\n")
  
  stopifnot(samp2 %in% norm_samp2)
  cat("-> all ", cond2, " in ", norm_cond2, ":\tok\n")
  
  if(cond1 == "wt") {
    stopifnot(setequal(c(samp1,samp2), norm_samp2))
    cat("-> setequal (", cond1, "+", cond2, ") "," and ", norm_cond2, ":\tok\n")
  }
}



