startTime <- Sys.time()

library(foreach)
library(doMC)
library(ggplot2)

options(scipen=100)

cat("> START: cmp_datasets_nbrTADs.R\n")
# Rscript cmp_datasets_nbrTADs.R

buildTable <- TRUE

SSHFS <- FALSE
# setDir <- ifelse(SSHFS, "/media/electron", "")
setDir <- ifelse(SSHFS, "~/media/electron", "")


# if(SSHFS) setwd("/media/electron/mnt/etemp/marie/Cancer_HiC_data_TAD_DA")
if(SSHFS) setwd("~/media/electron/mnt/etemp/marie/Cancer_HiC_data_TAD_DA")

source(file.path("utils_fct.R"))
source(file.path("datasets_settings.R"))

binSize <- 40000
binSizeKb <- binSize/1000

folderSuffix <- paste0("_", binSizeKb, "kb")
TADfilePattern <- "_final_domains.txt$"
TADfolder <- "FINAL_DOMAINS"

all_chromos <- c(paste0("chr", 1:23), "chrX")  # also used for ordering on the plots

xlabType <- "ds_label"
stopifnot(xlabType %in% c("ds1", "ds_label"))


### SELECT WHICH CELL LINES TO INCLUDE IN THE COMPARISONS
# in total (incl. pipCon): 30
cl_to_cmp <- c(
  "MCF-7",
  "ENCSR549MGQ_T47D",
  "MCF-7ENCSR549MGQ_T47D",
  
  "GSE105318_DLD1",                             
  
  "ENCSR079VIJ_G401",
  "ENCSR401TBQ_Caki2",
  "ENCSR079VIJ_G401ENCSR401TBQ_Caki2",
  
  "GSE105381_HepG2",                           
  
  "ENCSR444WCZ_A549",
  "NCI-H460",
  "ENCSR444WCZ_A549NCI-H460",
  
  "K562",
  
  "GSM2334834_U266_HindIII",
  "GSM2334832_RPMI-8226_HindIII",
  "GSM2334834_U266_HindIIIGSM2334832_RPMI-8226_HindIII",
  
  "ENCSR834DXR_SK-N-MC",
  "ENCSR105KFX_SK-N-DZ",
  # NB: no consensus because of metastatic sites
  
  "Panc1_rep12",
  
  "ENCSR346DCU_LNCaP",
  "GSE73782_PC3",             
  "ENCSR346DCU_LNCaPGSE73782_PC3",
  "GSE73782_PC3_ICE",             
  "ENCSR346DCU_LNCaPGSE73782_PC3_ICE",
  
  "ENCSR312KHQ_SK-MEL-5",
  "ENCSR862OGI_RPMI-7951",
  "ENCSR312KHQ_SK-MEL-5ENCSR862OGI_RPMI-7951",
  
  "GSE105194_spinal_cord",
  "GSE105194_cerebellum",   
  "GSE105194_spinal_cordGSE105194_cerebellum",
  
  "pipelineConsensus"
)
stopifnot(cl_to_cmp %in% cl_names)
stopifnot(dir.exists(paste0(cl_to_cmp, folderSuffix)))

registerDoMC(ifelse(SSHFS, 2, 40))

outFold <- file.path("CMP_DATASETS_NBRTADS")
dir.create(outFold, recursive=TRUE)

logFile <- file.path(outFold, "cmp_datasets_nbrTADs_logFile.txt")
if(!SSHFS) system(paste0("rm -f ", logFile))
if(SSHFS) logFile <- ""

plotType <- "svg"
widthBoxplot <- 10
heightBoxplot <- 7
myHeight <- ifelse(plotType == "png", 400, 7)
myWidth <- ifelse(plotType == "png", 400, 7)

cexPlot <- 1.2

txt <- paste0("!! Hard-coded buildTable:\t", as.character(buildTable), "\n")
printAndLog(txt, logFile)
txt <- paste0("!! Hard-coded bin size:\t", binSize, "\n")
printAndLog(txt, logFile)
txt <- paste0("!! Hard-coded # cell lines compared:\t", length(cl_to_cmp), "\n")
printAndLog(txt, logFile)
txt <- paste0("!! Hard-coded cell lines compared:\t", paste0(cl_to_cmp, collapse=","), "\n")
printAndLog(txt, logFile)

ds1=cl_to_cmp[1]
ds1="MCF-7ENCSR549MGQ_T47D"

if(buildTable){
all_nbr_dt <- foreach(ds1 = cl_to_cmp, .combine="rbind") %dopar% {
  
  folder1 <- paste0(ds1, folderSuffix)
  stopifnot(dir.exists(folder1))
  all_files1 <- list.files(file.path(folder1, TADfolder), pattern=TADfilePattern, recursive=FALSE, full.names = TRUE)
  stopifnot(length(all_files1) > 0)  
  stopifnot(length(all_files1) <= 23)

  cat(paste0("*** START: ", ds1,  "\n"))
    
  file1 = all_files1[1]
  chr_nbr_dt <- foreach(file1 = all_files1, .combine='rbind') %do% {
    
    chromo <- gsub(".+_(chr.+?)_.+", "\\1", basename(file1))
    stopifnot(length(chromo) == 1)
    stopifnot(chromo %in% all_chromos)
    
    cat(paste0("> ", ds1, "  - ", chromo, "\n"))

    dt1 <- read.delim(file1, stringsAsFactors = FALSE, header=F, col.names = c("chromo", "start", "end"))
    if(is.character(dt1[1,2]) & is.character(dt1[2,2])) {
      dt1 <- read.delim(file1, stringsAsFactors = FALSE, header=T, col.names = c("chromo", "start", "end"))
    }
    stopifnot(ncol(dt1) == 3)
    stopifnot(is.numeric(dt1[,2]))
    stopifnot(is.numeric(dt1[,3]))
    
    head(dt1, 2)
    
    chrEnd <- dt1[nrow(dt1), 3]
    chromoCover <- sum(dt1[,3] - dt1[,2] + 1)/chrEnd
    stopifnot(chromoCover >= 0 & chromoCover <= 1)
    
    data.frame(
      ds1 = ds1,
      chromo = chromo,
      chromoCover = chromoCover,
      nTADs = nrow(dt1),
      meanSizeTADs = mean(dt1$end - dt1$start + 1),
      stringsAsFactors = FALSE
    )
  }
  chr_nbr_dt
}

outFile <- file.path(outFold, "all_nbr_dt.Rdata")
save(all_nbr_dt, file = outFile)
cat(paste0("... written: ", outFile, "\n"))

} else {
  outFile <- file.path(outFold, "all_nbr_dt.Rdata")
  all_nbr_dt <- eval(parse(text = load(outFile)))
}


# take intersect chromos here
intersectChromos <- Reduce(intersect, by(all_nbr_dt$ds1, data=all_nbr_dt, FUN=function(x) unique(as.character(x$chromo))))

txt <- paste0("... # intersectChromos:\t", length(intersectChromos), "\n")
printAndLog(txt, logFile)
txt <- paste0("... intersectChromos:\t", paste0(intersectChromos, collapse=","), "\n")
printAndLog(txt, logFile)

nDS <- length(unique(all_nbr_dt$ds1))
all_nbr_dt$chromo <- as.character(all_nbr_dt$chromo)
all_nbr_dt <- all_nbr_dt[all_nbr_dt$chromo %in% intersectChromos, ]
stopifnot(nrow(all_nbr_dt) == (length(intersectChromos) * nDS))

all_nbr_dt$chromo <- factor(as.character(all_nbr_dt$chromo), levels = all_chromos)
stopifnot(!is.na(all_nbr_dt))

all_nbr_dt$ds_label <- sapply(as.character(all_nbr_dt$ds1), function(x) {   # !!! NEED THE AS.CHARACTER HERE !!!
  dslab <- as.character(cl_labs[x])
  stopifnot(length(dslab) == 1)
  if(is.na(dslab)) {
    paste0(names(cl_names[cl_names == x]))
  } else {
    paste0(names(cl_names[cl_names == x]), "\n(", dslab, ")")
  }
})

all_vars <- colnames(all_nbr_dt)[!colnames(all_nbr_dt) %in% c("ds1", "chromo", "ds_label")]

### VARIABLES TO PLOT
plot_tit <- c(
meanSizeTADs = "TAD mean size (bp)",
chromoCover = "Chromosome coverage",
nTADs = "Nbr of TADs"
)
stopifnot(all_vars %in% names(plot_tit))

var_to_plot = "nTADs"
for(var_to_plot in all_vars) {
  
  mytit <- plot_tit[var_to_plot]
  
  # TO PLOT THE ORIGINAL CL NAME
  mean_all_nbr_dt <- aggregate(as.formula(paste0(var_to_plot, " ~ ", xlabType)), FUN=mean, data = all_nbr_dt)
  mean_all_nbr_dt <- mean_all_nbr_dt[order(mean_all_nbr_dt[, var_to_plot], decreasing = TRUE),]
  all_nbr_dt[, xlabType] <- factor(as.character(all_nbr_dt[, xlabType]), levels = mean_all_nbr_dt[, xlabType])

  p_common <- ggplot(all_nbr_dt, aes_string(x = xlabType, y = var_to_plot)) + # plot cl ID
    geom_boxplot(outlier.shape=NA) +
    scale_x_discrete(name="")+
    scale_y_continuous(name=paste0(mytit),
                       breaks = scales::pretty_breaks(n = 10))+ 
    labs(colour  = "") +
#    ggtitle(label = paste0(var_to_plot))+
    ggtitle(label = paste0(mytit))+
    theme( # Increase size of axis lines
      # top, right, bottom and left
      # plot.margin = unit(c(1, 1, 4.5, 1), "lines"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size=10),
      panel.grid = element_blank(),
      axis.text.x = element_text( hjust=1,vjust = 0.5, size=12, angle = 90),
      axis.line.x = element_line(size = .2, color = "black"),
      axis.line.y = element_line(size = .3, color = "black"),
      axis.text.y = element_text(color="black", hjust=1,vjust = 0.5),
      axis.title.y = element_text(color="black", size=12),
      axis.title.x = element_text(color="black", size=12),
      panel.border = element_blank(),
      panel.background = element_rect(fill = "transparent"),
      legend.background =  element_rect(),
      legend.key = element_blank()
      ) 
  if(SSHFS) p_all
  p_dot <- p_common +  geom_jitter(aes(colour = chromo))
  if(SSHFS) p_dot
  p_txt <- p_common + geom_text(aes(label=chromo, colour=chromo, fontface="bold"),size=2.5, position = position_jitter(w = 0.3)) + guides(colour = "none")
  if(SSHFS) p_txt
  
  # outFile <- file.path(outFold, paste0(var_to_plot, "_boxplot_chromoDots.", plotType))
  # ggsave(plot=p_dot, file = outFile, width = widthBoxplot, height = heightBoxplot)
  # cat(paste0("... written: ", outFile, "\n"))
  
  outFile <- file.path(outFold, paste0(var_to_plot, "_boxplot_chromoLabs.", plotType))
  ggsave(plot=p_txt, file = outFile, width = widthBoxplot, height = heightBoxplot)
  cat(paste0("... written: ", outFile, "\n"))
}

##### NBR TADs AND TAD SIZE

myx <- all_nbr_dt$nTADs
myy <- log10(all_nbr_dt$meanSizeTADs)

myxlab <- "# of TADs"
myylab <- "mean TAD size [log10]"

myTit <- "TAD size vs. nbr of TADs"

mySub <- paste0("(# DS = ", nDS, ")")

outFile <- file.path(outFold, paste0("nbrTADs_vs_TADsize_scatterplot.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))

densplot(x=myx,
         y=myy,
         xlab=myxlab,
         ylab=myylab,
         pch = 16, cex = 0.7,
         cex.lab = cexPlot, cex.axis = cexPlot,
         main = myTit
)
mtext(side=3, text = mySub)
add_curv_fit(x = myx,
             y = myy,
             withR2 = FALSE, lty=2, col="darkgray")

addCorr(x=myx,
        y=myy,
        bty="n",
        legPos="bottomright")
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))




######################################################################################
cat(paste0("... written: ", logFile, "\n"))
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))






