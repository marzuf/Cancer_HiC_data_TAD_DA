# COMPARE TO THE SAME FILE IN EZH2_final_MAPQ, HERE TAD HEADER FILE = TRUE


options(scipen=100)

startTime <- Sys.time()

suppressPackageStartupMessages(library(optparse, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

option_list = list(
  
  make_option(c("-f", "--feature_file"), type="character", default=NULL,
              help="input feature (gene) file", metavar="character"),            
  
  make_option(c("-m", "--matrix_file"), type="character", default=NULL,
              help="input matrix file", metavar="character"),              
  
  make_option(c("-s", "--start_matrix"), type="integer", default=NULL,
              help="draw from start_matrix (in bp !)", metavar="character"),              
  
  make_option(c("-e", "--end_matrix"), type="integer", default=NULL,
              help="draw to end_matrix (in bp !)", metavar="character"),              
  
  make_option(c("-o", "--output_file"), type="character", default=NULL,
              help="path to output file", metavar="character"),              
  
  make_option(c("-k", "--col_to_skip"), type="integer", default=NULL,
              help="columns to skip", metavar="integer"),
  
  make_option(c("-c", "--chromo"), type="character", default=NULL,
              help="chromosome to draw", metavar="character"),
  
  make_option(c("-b", "--bin_size"), type="integer", default=NULL,
              help="binning size", metavar="integer")
  
);

opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);

if(is.null(opt$matrix_file) |  is.null(opt$bin_size) | 
   is.null(opt$output_file) ) {
  stop("Missing arguments \n")
}

chromo <- opt$chromo


featureFile <- opt$feature_file

matrixFile <- opt$matrix_file
binSize <- opt$bin_size
skipcol <- ifelse(is.null(opt$col_to_skip), 3, opt$col_to_skip)
start_matrix <- opt$start_matrix
end_matrix <- opt$end_matrix

outFile <- opt$output_file
system(paste0("mkdir -p ", dirname(outFile)))

stopifnot(file.exists(matrixFile))
if(!is.null(featureFile)) stopifnot(file.exists(featureFile))

system(paste0("mkdir -p ", dirname(outFile)))

plotType <- gsub(".+\\.(.+?)$", "\\1", basename(outFile))
myHeight <- ifelse(plotType == "pdf" | plotType == "svg", 7, 480)
myWidth <- myHeight

imageColPalette <- colorRampPalette(c("blue", "red"))( 12 )

matrixFormat <- "domaincaller"

########################################## HARD-CODED PARAMETERS
matrixHeader <- FALSE
featureHeader <- FALSE

featureCol <- "cyan"


#### DROP THE FIRST COLUMNS OF THE MATRIX
cat(paste0("... load matrix data\t", Sys.time(), "\t"))

if(matrixFormat == "dekker") {
    matrixDT <- read.delim(matrixFile, header=T, skip = 1, check.names = F)
    cat(paste0(Sys.time(), "\n"))
    rownames(matrixDT) <- matrixDT[,1]
    matrixDT[,1] <- NULL
    stopifnot(ncol(matrixDT) == nrow(matrixDT) + skipcol)
    stopifnot(!is.na(colnames(matrixDT)))
    stopifnot(!is.na(rownames(matrixDT)))
    if(skipcol > 0)
      matrixDT <- matrixDT[,-c(1:skipcol)]
    stopifnot(colnames(matrixDT) == rownames(matrixDT))
    stopifnot(nrow(matrixDT) == ncol(matrixDT) )
    rownames(matrixDT) <- colnames(matrixDT) <- NULL
} else {
    matrixDT <- read.delim(matrixFile, header=matrixHeader, stringsAsFactors = FALSE)
    cat(paste0(Sys.time(), "\n"))
    stopifnot(ncol(matrixDT) == nrow(matrixDT) + skipcol)
    if(skipcol > 0)
      matrixDT <- matrixDT[,-c(1:skipcol)]
    stopifnot(nrow(matrixDT) == ncol(matrixDT) )
}


cat("... discard data don't want to plot\n")

#### PREPARE THE MATRIX - SELECT FROM THE MATRIX THE AREA WE WANT TO PLOT
if(is.null(start_matrix)) {
  start_matrix <- 1
} else {
  # convert the start limit in bp to bin
  start_matrix <- floor(start_matrix/binSize) + 1
}

if(start_matrix > nrow(matrixDT)) {
  stop("... want to start plotting after the end of the matrix!\n")
}

if(is.null(end_matrix)) {
  end_matrix <- nrow(matrixDT)
} else {
  end_matrix <-  ceiling(end_matrix/binSize) 
  if(end_matrix > ncol(matrixDT)){
    cat("! WARNING: wanted end position is after end of the data, will plot up to the end\n")
    end_matrix <- ncol(matrixDT)
  }
}
stopifnot(end_matrix >= start_matrix)
stopifnot(start_matrix > 0 & end_matrix <= ncol(matrixDT))
cat("... will draw from bin:\t", start_matrix, "\tto:\t", end_matrix , "(inclusive)\n")

matrixDT <- matrixDT[start_matrix:end_matrix, start_matrix:end_matrix]
# revert the matrix to have the plot from topleft to bottom right
drawMatrixDT <- t(matrixDT)[,nrow(matrixDT):1]

#### PREPARE THE TADs - ADJUST POSITIONS
shift_bin <- start_matrix - 1

do.call(plotType, list(outFile, height=myHeight, width=myWidth))
totBin <- nrow(matrixDT) + 1
axLab <- seq(1.5, length.out=nrow(matrixDT))
# image(x=axLab, y=axLab, as.matrix(drawMatrixDT),
#       xlab="", ylab="",
#       xaxt = "n", yaxt="n")
cat("... draw the image\n")
image(x=axLab, y=axLab, as.matrix(log10(drawMatrixDT+0.001)),
      xlab="", ylab="",
      xaxt = "n", yaxt="n",
      col = imageColPalette)

title(paste0(chromo, " - ", 1+binSize*(start_matrix-1), "(", start_matrix, "):", end_matrix*binSize, "(", end_matrix,")"))

### add starts for the genes if provided
if(!is.null(featureFile)){
  cat("... add feature segments \n")
  featureDT <- read.delim(featureFile, header=featureHeader, stringsAsFactors = FALSE)
  if(ncol(featureDT) == 3){
    colnames(featureDT) <- c("chromo", "start", "end")
    labelFeature <- FALSE
  } else if(ncol(featureDT) == 4){
    colnames(featureDT) <- c("chromo", "start", "end", "gene")
  } else{
    stop("unknown format feature file\n")
  }
  
  featureDT <- featureDT[featureDT$chromo == chromo,]
  
  if(nrow(featureDT) > 0){
    for(i in 1:nrow(featureDT)) {
      firstBin <- floor(featureDT$start[i]/binSize)+1 -shift_bin
      lastBin <- ceiling(featureDT$end[i]/binSize) -shift_bin
      stopifnot(lastBin >= firstBin)
      for(feat_bin in firstBin:lastBin){
        my_xpos <- (feat_bin + (feat_bin+1))*0.5
        my_ypos <-  (totBin-feat_bin+1 + totBin -feat_bin)*0.5
        points(x=my_xpos, y=my_ypos, pch=16, cex = 1, adj=0.5, col = featureCol)
      }
    }
  }
}
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))

