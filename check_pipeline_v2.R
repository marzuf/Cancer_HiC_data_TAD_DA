setwd("~/media/electron/mnt/etemp/marie/Cancer_HiC_data_TAD_DA")

plotType <- "svg"
myHeightDensity <- ifelse(plotType == "png", 400, 7)
myWidthDensity <- ifelse(plotType == "png", 600, 10)
myHeight <- myWidth <- myHeightDensity

plotCex <- 1.2

outFold <- "CHECK_PIPELINE_V2"

step_2v2 <- TRUE

step_5v2 <- TRUE

pipOutFolder <- "PIPELINE/OUTPUT_FOLDER"

hicds <- "ENCSR079VIJ_G401_40kb"
exprds <- "TCGAkich_norm_kich"

mainFolder <- file.path(pipOutFolder, hicds, exprds)
stopifnot(dir.exists(mainFolder))

################################################
################################################ STEP 2v2
################################################
script2v2_name <- "2v2_runWilcoxonTAD"
script3_name <- "3_runMeanTADLogFC"

signifThresh1 <- 0.05
signifThresh2 <- 0.01

#================> density plots (Wilcox test pval)
mytit <- paste0(hicds, " - ", exprds)
mysub <- paste0("FPKM gene avgd. TAD paired Wilcox. test")
myxlab <- paste0("adj. wilcox. p-val.")

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

wilcox_pvals <- unlist(lapply(wilcox_pairedTAD_meanExpr_fpkm, function(x) x[["wilcoxTest_pval"]]))
adj_wilcox_pvals <- p.adjust(wilcox_pvals, method="BH")

plot(density(adj_wilcox_pvals), 
     main = paste0(mytit),
     xlab = myxlab,
     cex.lab = plotCex,
     cex.axis = plotCex)
mtext(text = paste0(mysub), side=3)

nSignif1 <- sum(adj_wilcox_pvals <= signifThresh1)
nSignif2 <- sum(adj_wilcox_pvals <= signifThresh2)

legend("topleft", 
       legend = c(paste0("# TADs = ", length(adj_wilcox_pvals)),
                  paste0("# <=", signifThresh1, "=", nSignif1),
                  paste0("# <=", signifThresh2, "=", nSignif2)),
       bty="n"
       )

#================================
#================> density plots (Wilcox test pval)
#================================
mytit <- paste0(hicds, " - ", exprds)
mysub <- paste0("FPKM gene avgd. TAD paired Wilcox. test")
myxlab <- paste0("Wilcox. W stat.")

wilcox_wstat <- unlist(lapply(wilcox_pairedTAD_meanExpr_fpkm, function(x) x[["wilcoxTest_stat"]]))

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


#================================
#================> comparison with v1
#================================
v1_meanLogFC_file <- file.path(mainFolder, script3_name, "all_meanLogFC_TAD.Rdata")
load(v1_meanLogFC_file)
head(all_meanLogFC_TAD)
stopifnot(setequal(names(all_meanLogFC_TAD),  names(adj_wilcox_pvals)))

tad_list <- names(all_meanLogFC_TAD)
stopifnot(!duplicated(tad_list))

mytit <- paste0(hicds, " - ", exprds)
mysub <- paste0("Wilcox. test adj. pval vs. meanLogFC")
myylab <- paste0("Wilcox. test adj. pval (v2)")
myxlab <- paste0("TAD mean logFC (v1)")

# densplot(y = adj_wilcox_pvals[tad_list],
#          x = all_meanLogFC_TAD[tad_list],
#      main = paste0(mytit),
#      xlab = myxlab,
#      ylab = myylab,
#      cex.lab = plotCex,
#      cex.axis = plotCex)
# mtext(text = paste0(mysub), side=3)
# 
# legend("topleft", 
#        legend = c(paste0("# TADs = ", length(adj_wilcox_pvals))),
#        bty="n"
# )
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


################################################
################################################ STEP 5v2
################################################

script5v2_name <- "5v2_runPermutations"

#================================
#================> density plot: number of exactly same TAD during permutations
#================================

v5_nbrSameTAD_file <- file.path(mainFolder, script5v2_name, "nbrSameTAD_obs_permut.Rdata")
load(v5_nbrSameTAD_file)
head(nbrSameTAD_obs_permut)

stopifnot(length(unique(unlist(lapply(nbrSameTAD_obs_permut, length)))) == 1)

stopifnot(length(unique(unlist(lapply(nbrSameTAD_obs_permut, length)))) == 1)

n1 <- names(nbrSameTAD_obs_permut[[1]])
stopifnot(unlist(lapply(nbrSameTAD_obs_permut, function(x) setequal(names(x), n1))))
stopifnot(!duplicated(n1))
nTotTADs <- length(n1)

mytit <- paste0(hicds, " - ", exprds)
mysub <- paste0("Ratio same TADs during permut.")
myxlab <- paste0("ratio exactly same TAD")

nSameTADs <- unlist(lapply(nbrSameTAD_obs_permut, sum))

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





