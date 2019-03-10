SSHFS <- TRUE
setDir <- ifelse(SSHFS, "~/media/electron", "")
setwd(setDir)

plotType <- "svg"
myHeightDensity <- ifelse(plotType == "png", 400, 7)
myWidthDensity <- ifelse(plotType == "png", 600, 10)
myHeight <- myWidth <- myHeightDensity
myWidthGG <- myWidthDensity
plotCex <- 1.2

outFold <- file.path("CHECK_PIPELINE_V2_ALLDS")
dir.create(outFold, recursive=TRUE)

step_2v2 <- TRUE
step_5v2 <- TRUE
step_7v2 <- TRUE

pipOutFolder <- "PIPELINE/OUTPUT_FOLDER"

all_hicexpr_ds <- unname(unlist(sapply(list.files(pipOutFolder, full.names = TRUE), function(x) file.path(basename(x),  list.files(x)))))
stopifnot(dir.exists(file.path(pipOutFolder, all_hicexpr_ds)))

build_cmp_v2_allDS_data=TRUE

ds=all_hicexpr_ds[1]

hicds <- "ENCSR079VIJ_G401_40kb"
exprds <- "TCGAkich_norm_kich"

mainFolder <- file.path(pipOutFolder, hicds, exprds)
stopifnot(dir.exists(mainFolder))

source("utils_fct.R")



script2v2_name <- "2v2_runWilcoxonTAD"
script3_name <- "3_runMeanTADLogFC"


script5v2_name <- "5v2_runPermutations" # 44

script5v2_name <- "5v2_runPermutations" # 44


script7_name <- "7_runPermutationsMeanTADCorr"
script7v2_name <- "7v2_runPermutationsMeanTADCorr"
script4_name <- "4_runMeanTADCorr"

script10_name <- "10_runEmpPvalMeanTADCorr"
script10v2_name <- "10v2_runEmpPvalMeanTADCorr"



#============================= <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
#============================= <<<<<<<<<<<<<<<<<<<<<<<<<<<<< COLLECT DATA
#============================= <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 



if(build_cmp_v2_allDS_data){
  
  cat("... start building signifTADs_allDS_data \n")
  cmp_v2_allDS_data <- foreach(ds = all_hicexpr_ds) %dopar% {
    
    cat("... start: ", ds, "\n")
    
    hicds <- dirname(ds)
    exprds <- basename(ds)
    
    ########################################################################################################
    ########################################################################################################
    ########################################################################################################
    
    
    
    
    ################################################
    ################################################ STEP 2v2 => wilcox stat obs. data
    ################################################

    
    
    wilcoxTestFile <- file.path(mainFolder, script2v2_name, "wilcox_pairedTAD_meanExpr_wilcoxStat.Rdata")
    stopifnot(file.exists(wilcoxTestFile))
    load(wilcoxTestFile)
    head(wilcox_pairedTAD_meanExpr_wilcoxStat)
    head(wilcox_pairedTAD_meanExpr_wilcoxStat)
    # chr10_TAD1  chr10_TAD2  chr10_TAD8  chr10_TAD9 chr10_TAD10 chr10_TAD11 
    # 10          21          20          10          15          10 
    
    wilcoxStatFile <- file.path(mainFolder, script2v2_name, "wilcox_pairedTAD_meanExpr_fpkm.Rdata")
    stopifnot(file.exists(wilcoxStatFile))
    load(wilcoxStatFile)
    head(wilcox_pairedTAD_meanExpr_fpkm)
    # $chr10_TAD11
    # $chr10_TAD11$wilcoxTest_pval
    # [1] 0.125
    # 
    # $chr10_TAD11$wilcoxTest_stat
    # [1] 10
    
    wilcox_wstat <- unlist(lapply(wilcox_pairedTAD_meanExpr_fpkm, function(x) x[["wilcoxTest_stat"]]))
    
    
    wilcox_pvals <- unlist(lapply(wilcox_pairedTAD_meanExpr_fpkm, function(x) x[["wilcoxTest_pval"]]))
    adj_wilcox_pvals <- p.adjust(wilcox_pvals, method="BH")
    
    nSignif1 <- sum(adj_wilcox_pvals <= signifThresh1)
    nSignif2 <- sum(adj_wilcox_pvals <= signifThresh2)
    

        
    #================================
    #================> density plots (Wilcox test pval)
    #================================
    wilcox_wstat <- unlist(lapply(wilcox_pairedTAD_meanExpr_fpkm, function(x) x[["wilcoxTest_stat"]]))
    

    #================================
    #================> comparison with v1 
    #================================
    v1_meanLogFC_file <- file.path(mainFolder, script3_name, "all_meanLogFC_TAD.Rdata")
    load(v1_meanLogFC_file)
    head(all_meanLogFC_TAD)
    stopifnot(setequal(names(all_meanLogFC_TAD),  names(adj_wilcox_pvals)))
    
    tad_list <- names(all_meanLogFC_TAD)
    stopifnot(!duplicated(tad_list))
    

    ################################################
    ################################################ STEP 5v2 => permutations for intraCorr
    ################################################
    
    
    #================================
    #================> density plot: number of exactly same TAD during permutations
    #================================
    
    v2_nbrSameTAD_file <- file.path(mainFolder, script5v2_name, "nbrSameTAD_obs_permut.Rdata")
    load(v2_nbrSameTAD_file)
    head(nbrSameTAD_obs_permut)
    
    stopifnot(length(unique(unlist(lapply(nbrSameTAD_obs_permut, length)))) == 1)
    
    stopifnot(length(unique(unlist(lapply(nbrSameTAD_obs_permut, length)))) == 1)
    
    n1 <- names(nbrSameTAD_obs_permut[[1]])
    stopifnot(unlist(lapply(nbrSameTAD_obs_permut, function(x) setequal(names(x), n1))))
    stopifnot(!duplicated(n1))
    nTotTADs <- length(n1)
    

    nSameTADs <- unlist(lapply(nbrSameTAD_obs_permut, sum))
    

    ################################################
    ################################################ STEP 7v2 => intraCorr for permutDT
    ################################################
    
    v2_permutMeanTADcorr_file <- file.path(mainFolder, script7v2_name, "meanCorr_permDT.Rdata")
    load(v2_permutMeanTADcorr_file)
    meanCorr_permDT[1:5,1:5]
    meanCorr_permDT_v2 <- meanCorr_permDT
    
    
    v1_permutMeanTADcorr_file <- file.path(mainFolder, script7_name, "meanCorr_permDT.Rdata")
    load(v1_permutMeanTADcorr_file)
    meanCorr_permDT[1:5,1:5]
    meanCorr_permDT_v1 <- meanCorr_permDT
    
    # in the mean time:
    meanCorr_permDT_v2 <- eval(parse(text = load("luad_kras_egfr_foo_7v2_meanCorr_permDT.Rdata")))
    meanCorr_permDT_v1 <- eval(parse(text = load("luad_kras_egfr_foo_7_meanCorr_permDT.Rdata")))
    
    stopifnot(dim(meanCorr_permDT_v2) == dim(meanCorr_permDT) )
    
    stopifnot(rownames(meanCorr_permDT_v1) == rownames(meanCorr_permDT_v2))
    
    stopifnot(!is.na(meanCorr_permDT_v1))
    stopifnot(!is.na(meanCorr_permDT_v2))
    
    #================================
    #================> densplot comparison with v1 
    #================================

    tad_avg_v1 <- rowMeans(meanCorr_permDT_v1)
    tad_avg_v2 <- rowMeans(meanCorr_permDT_v2)
    stopifnot(names(tad_avg_v1) == names(tad_avg_v2) )
    

    #================================
    #================> multi density plot:  comparison with v1 
    #================================
    

    #================================
    #================> multi density plot:  comparison with v1 and observed 
    #================================
    obs_corr_file <-  file.path(mainFolder, script4_name, "all_meanCorr_TAD.Rdata")
    obs_corr <- eval(parse(text = load(obs_corr_file)))
    

    # > x="PIPELINE/OUTPUT_FOLDER/GSE105381_HepG2_40kb/TCGAlihc_norm_lihc/4_runMeanTADCorr/all_meanCorr_TAD.Rdata""
    
    ################################################
    ################################################ STEP 10v2 => emp. pval intraCorr 
    ################################################
    # > emp_pval_meanCorr_v1_file ="PIPELINE/OUTPUT_FOLDER/ENCSR079VIJ_G401_40kb/TCGAkich_norm_kich/10_runEmpPvalMeanTADCorr/emp_pval_meanCorr.Rdata"
    # > emp_pval_meanCorr_v2_file ="PIPELINE/OUTPUT_FOLDER/ENCSR079VIJ_G401_40kb/TCGAkich_norm_kich/10v2_runEmpPvalMeanTADCorr/emp_pval_meanCorr.Rdata"
    
    
    emp_pval_meanCorr_v1_file <- file.path(mainFolder, script10_name, "emp_pval_meanCorr.Rdata")
    emp_pval_meanCorr_v2_file <- file.path(mainFolder, script10v2_name, "emp_pval_meanCorr.Rdata")
    
    emp_pval_meanCorr_v1 <- eval(parse(text = load(emp_pval_meanCorr_v1_file)))
    emp_pval_meanCorr_v2 <- eval(parse(text = load(emp_pval_meanCorr_v2_file)))
    
    adj_emp_pval_meanCorr_v1 <- p.adjust(emp_pval_meanCorr_v1, method="BH")
    adj_emp_pval_meanCorr_v2 <- p.adjust(emp_pval_meanCorr_v2, method="BH")
    
    stopifnot(length(adj_emp_pval_meanCorr_v1) == length(adj_emp_pval_meanCorr_v2) )
    stopifnot(setequal(names(adj_emp_pval_meanCorr_v1), names(adj_emp_pval_meanCorr_v2)))
    stopifnot( names(adj_emp_pval_meanCorr_v1) == names(adj_emp_pval_meanCorr_v2) )
    
    

    #================================
    #================> densplot comparison with v1 
    #================================
    

    #================================
    #================> densplot comparison with v1 - log10
    #================================
    

    #================================
    #================> multi density plot:  comparison with v1 
    #================================
    
    #================================
    #================> signif by threshold
    #================================
    
    
    pval_step <- 0.01
    
    pval_signif_thresh_seq <- seq(from=0, to=1, by=pval_step) 
    
    ratioSignif_by_thresh_v1 <- sapply(pval_signif_thresh_seq, function(p_thresh) {
      mean(adj_emp_pval_meanCorr_v1 <= p_thresh) 
    })
    names(ratioSignif_by_thresh_v1) <- pval_signif_thresh_seq
    head(ratioSignif_by_thresh_v1)
    
    ratioSignif_by_thresh_v2 <- sapply(pval_signif_thresh_seq, function(p_thresh) {
      mean(adj_emp_pval_meanCorr_v2 <= p_thresh) 
    })
    names(ratioSignif_by_thresh_v2) <- pval_signif_thresh_seq
    head(ratioSignif_by_thresh_v2)
  

    
    outFile <- file.path(outFold, "cmp_v2_allDS_data.Rdata")
    save(cmp_v2_allDS_data, file = outFile)
    cat(paste0("... written: ", outFile, "\n"))
    

  } # end foreach all_hicexpr_ds
} else {
  
  outFile <- file.path(outFold, "cmp_v2_allDS_data.Rdata")
  cmp_v2_allDS_data <- eval(parse(text=load(outFile)))
}


#============================= <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
#============================= <<<<<<<<<<<<<<<<<<<<<<<<<<<<< VARIABLES TO RETRIEVE
#============================= <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 

################################################
################################################ STEP 2v2 => wilcox stat obs. data
################################################
adj_wilcox_pvals
#================================
#================> density plots (Wilcox test pval)
#================================
wilcox_wstat
#================================
#================> comparison with v1 
#================================
(y = -log10(adj_wilcox_pvals[tad_list]),
 x = all_meanLogFC_TAD[tad_list],
 
 ################################################
 ################################################ STEP 5v2 => permutations for intraCorr
 ################################################
 
 
 #================================
 #================> density plot: number of exactly same TAD during permutations
 #================================
 nSameTADs/nTotTADs
 
 ################################################
 ################################################ STEP 7v2 => intraCorr for permutDT
 ################################################
 
 
 
 #================================
 #================> densplot comparison with v1 
 #================================
 meanCorr_permDT_v1
 meanCorr_permDT_v2
 
 
 tad_avg_v1 
 
 tad_avg_v2 
 
 
 #================================
 #================> multi density plot:  comparison with v1 
 #================================
 
 meanCorr_permDT_v1
 meanCorr_permDT_v2
 
 
 #================================
 #================> multi density plot:  comparison with v1 and observed 
 #================================
 
 meanCorr_permDT_v1
 meanCorr_permDT_v2
 
 # > x="PIPELINE/OUTPUT_FOLDER/GSE105381_HepG2_40kb/TCGAlihc_norm_lihc/4_runMeanTADCorr/all_meanCorr_TAD.Rdata""
 
 ################################################
 ################################################ STEP 10v2 => emp. pval intraCorr 
 ################################################
 
 #================================
 #================> densplot comparison with v1 
 #================================
 x = adj_emp_pval_meanCorr_v1,
 y = adj_emp_pval_meanCorr_v2,
 
 
 #================================
 #================> densplot comparison with v1 - log10
 #================================
 
 
 log10(adj_emp_pval_meanCorr_v1),
 y = log10(adj_emp_pval_meanCorr_v2),
 
 
 
 
 #================================
 #================> multi density plot:  comparison with v1 
 #================================
 v1=as.numeric(unlist(adj_emp_pval_meanCorr_v1)),
 v2=as.numeric(unlist(adj_emp_pval_meanCorr_v2))
),


#================================
#================> signif by threshold
#================================

pval_step <- 0.01

pval_signif_thresh_seq <- seq(from=0, to=1, by=pval_step) 

xlim=range(c(ratioSignif_by_thresh_v1, ratioSignif_by_thresh_v2)),
ylim=range(c(ratioSignif_by_thresh_v1, ratioSignif_by_thresh_v2)), 

lines(x = pval_signif_thresh_seq, y = ratioSignif_by_thresh_v1, col="black")
lines(x = pval_signif_thresh_seq, y = ratioSignif_by_thresh_v2, col="red")




#============================= <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
#============================= <<<<<<<<<<<<<<<<<<<<<<<<<<<<< PLOTS
#============================= <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 






################################################
################################################ STEP 2v2 => wilcox stat obs. data
################################################




outFile <- file.path(outFold, paste0("2v2_adj_wilcox_pvals_density.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot(density(adj_wilcox_pvals), 
     main = paste0(mytit),
     xlab = myxlab,
     cex.lab = plotCex,
     cex.axis = plotCex)
mtext(text = paste0(mysub), side=3)
legend("topleft", 
       legend = c(paste0("# TADs = ", length(adj_wilcox_pvals)),
                  paste0("# <=", signifThresh1, "=", nSignif1),
                  paste0("# <=", signifThresh2, "=", nSignif2)),
       bty="n"
)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

#================================
#================> density plots (Wilcox test pval)
#================================
mytit <- paste0(hicds, " - ", exprds)
mysub <- paste0("FPKM gene avgd. TAD paired Wilcox. test")
myxlab <- paste0("Wilcox. W stat.")


outFile <- file.path(outFold, paste0("2v2_adj_wilcox_wstat_density.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot(density(wilcox_wstat), 
     main = paste0(mytit),
     xlab = myxlab,
     cex.lab = plotCex,
     cex.axis = plotCex)
mtext(text = paste0(mysub), side=3)

legend("topleft", 
       legend = c(paste0("# TADs = ", length(adj_wilcox_pvals))),
       bty="n"
)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

#================================
#================> comparison with v1 
#================================

mytit <- paste0(hicds, " - ", exprds)
mysub <- paste0("Wilcox. test adj. pval vs. meanLogFC")
myylab <- paste0("Wilcox. test adj. pval (v2)")
myxlab <- paste0("TAD mean logFC (v1)")


outFile <- file.path(outFold, paste0("2v2_adj_wilcox_pvals_densplot_log10.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidthGG))
densplot(y = -log10(adj_wilcox_pvals[tad_list]),
         x = all_meanLogFC_TAD[tad_list],
         main = paste0(mytit),
         xlab = myxlab,
         ylab = paste0(myylab, " [-log10]"),
         cex.lab = plotCex,
         cex.axis = plotCex)
mtext(text = paste0(mysub), side=3)

legend("topleft", 
       legend = c(paste0("# TADs = ", length(adj_wilcox_pvals))),
       bty="n"
)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


################################################
################################################ STEP 5v2 => permutations for intraCorr
################################################


#================================
#================> density plot: number of exactly same TAD during permutations
#================================


mytit <- paste0(hicds, " - ", exprds)
mysub <- paste0("Ratio same TADs during permut.")
myxlab <- paste0("ratio exactly same TAD")


outFile <- file.path(outFold, paste0("5v2_ratio_sameTADs_density.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot(density(nSameTADs/nTotTADs), 
     main = paste0(mytit),
     xlab = myxlab,
     cex.lab = plotCex,
     cex.axis = plotCex)
mtext(text = paste0(mysub), side=3)

legend("topleft", 
       legend = c(paste0("# TADs = ", length(n1)),
                  paste0("mean same TAD = ", mean(nSameTADs))),
       bty="n"
)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

################################################
################################################ STEP 7v2 => intraCorr for permutDT
################################################



#================================
#================> densplot comparison with v1 
#================================
mytit <- paste0(hicds, " - ", exprds)
myxlab <- "permut. v1"
myylab <- "permut. v2"
mysub <- "permut. avg. meanCorr"

tad_avg_v1 <- rowMeans(meanCorr_permDT_v1)
tad_avg_v2 <- rowMeans(meanCorr_permDT_v2)
stopifnot(names(tad_avg_v1) == names(tad_avg_v2) )

outFile <- file.path(outFold, paste0("7v2_permutIntraCorr_densplot.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidthGG))
densplot(
  x = tad_avg_v1,
  y = tad_avg_v2,
  main = paste0(mytit),
  xlab = myxlab,
  ylab = myylab,
  cex.lab = plotCex,
  cex.axis = plotCex)
mtext(text = paste0(mysub), side=3)

legend("topleft", 
       legend = c(paste0("# TADs = ", length(tad_avg_v1))),
       bty="n"
)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

#================================
#================> multi density plot:  comparison with v1 
#================================

outFile <- file.path(outFold, paste0("7v2_permutIntraCorr_multiDens.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidthGG))

plot_multiDens(size_list = list(
  v1=as.numeric(unlist(meanCorr_permDT_v1)),
  v2=as.numeric(unlist(meanCorr_permDT_v2))
),
plotTit = mytit,
my_xlab = "permut. meanCorr")
mtext(text = paste0(mysub), side=3)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

#================================
#================> multi density plot:  comparison with v1 and observed 
#================================

outFile <- file.path(outFold, paste0("7v2_permutIntraCorr_with_obs_multiDens.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidthGG))

plot_multiDens(size_list = list(
  v1=as.numeric(unlist(meanCorr_permDT_v1)),
  v2=as.numeric(unlist(meanCorr_permDT_v2)),
  obs=obs_corr
),
plotTit = mytit,
my_xlab = "permut. meanCorr")
mtext(text = paste0(mysub), side=3)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

# > x="PIPELINE/OUTPUT_FOLDER/GSE105381_HepG2_40kb/TCGAlihc_norm_lihc/4_runMeanTADCorr/all_meanCorr_TAD.Rdata""

################################################
################################################ STEP 10v2 => emp. pval intraCorr 
################################################
# > emp_pval_meanCorr_v1_file ="PIPELINE/OUTPUT_FOLDER/ENCSR079VIJ_G401_40kb/TCGAkich_norm_kich/10_runEmpPvalMeanTADCorr/emp_pval_meanCorr.Rdata"
# > emp_pval_meanCorr_v2_file ="PIPELINE/OUTPUT_FOLDER/ENCSR079VIJ_G401_40kb/TCGAkich_norm_kich/10v2_runEmpPvalMeanTADCorr/emp_pval_meanCorr.Rdata"

mytit <- paste0(hicds, " - ", exprds)
myxlab <- "adj. emp. pval. meanCorr v1"
myylab <- "adj. emp. pval. meanCorr v2"
mysub <- "adj. emp. pval. meanCorr"

#================================
#================> densplot comparison with v1 
#================================

outFile <- file.path(outFold, paste0("10v2_emp_pval_intraCorr_densplot.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidthGG))

densplot(
  x = adj_emp_pval_meanCorr_v1,
  y = adj_emp_pval_meanCorr_v2,
  main = paste0(mytit),
  xlab = myxlab,
  ylab = myylab,
  cex.lab = plotCex,
  cex.axis = plotCex)
mtext(text = paste0(mysub), side=3)

legend("topleft", 
       legend = c(paste0("# TADs = ", length(tad_avg_v1))),
       bty="n"
)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

#================================
#================> densplot comparison with v1 - log10
#================================

outFile <- file.path(outFold, paste0("10v2_emp_pval_intraCorr_densplot_log10.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidthGG))

densplot(
  x = log10(adj_emp_pval_meanCorr_v1),
  y = log10(adj_emp_pval_meanCorr_v2),
  main = paste0(mytit),
  xlab = paste0(myxlab, " [log10]"),
  ylab = paste0(myylab, " [log10]"),
  cex.lab = plotCex,
  cex.axis = plotCex)
mtext(text = paste0(mysub), side=3)

legend("topleft", 
       legend = c(paste0("# TADs = ", length(tad_avg_v1))),
       bty="n"
)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

#================================
#================> multi density plot:  comparison with v1 
#================================

outFile <- file.path(outFold, paste0("10v2_emp_pval_intraCorr_multidens.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidthGG))

plot_multiDens(size_list = list(
  v1=as.numeric(unlist(adj_emp_pval_meanCorr_v1)),
  v2=as.numeric(unlist(adj_emp_pval_meanCorr_v2))
),
plotTit = mytit,
my_xlab = "adj. emp. p-val meanCorr")
mtext(text = paste0(mysub), side=3)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

#================================
#================> signif by threshold
#================================

pval_step <- 0.01

pval_signif_thresh_seq <- seq(from=0, to=1, by=pval_step) 

myTit <- paste0("ratio signif. TADs vs. p-val. thresh.")
mySub <- paste0("(pval_step = ", pval_step, ")")
myylab <- paste0("ratio signif. TADs")
myxlab <- paste0("adj. emp. p-val thresh")

outFile <- file.path(outFold, paste0("10v2_emp_pval_intraCorr_nSignif_by_thresh.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))

plot(NULL,
     xlim=range(c(ratioSignif_by_thresh_v1, ratioSignif_by_thresh_v2)),
     ylim=range(c(ratioSignif_by_thresh_v1, ratioSignif_by_thresh_v2)), 
     cex.lab = plotCex, cex.axis = plotCex,
     xlab = myxlab,
     main = myTit)
mtext(text = mySub, side = 3)
lines(x = pval_signif_thresh_seq, y = ratioSignif_by_thresh_v1, col="black")
lines(x = pval_signif_thresh_seq, y = ratioSignif_by_thresh_v2, col="red")
legend("topleft", 
       lty=c(-1,1,1),
       col = c("black", "red"),
       legend = c(paste0("# TADs=", length(adj_emp_pval_meanCorr_v1)), "v1", "v2"),
       bty="n")
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))




