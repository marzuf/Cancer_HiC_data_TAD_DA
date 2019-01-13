#!/usr/bin/env Rscript

# just update because with the other versions the input files are not 3 columns, have also the idx and score columns as well as an header

options(scipen=100)

library(foreach)
library(optparse,verbose=F, quietly=T)

option_list = list(
  make_option(c("-f", "--gene_file"), type="character", default=NULL,
              help="path to file with genes", metavar="character"),
  make_option(c("-t", "--tad_file"), type="character", default=NULL,
              help="path to file with TAD", metavar="character"),
  make_option(c("-b", "--bin_size"), type="integer", default=NULL,
              help="bin size", metavar="character"),
  make_option(c("-c", "--chromo"), type="character", default=NULL,
              help="chromo name", metavar="character"),
  make_option(c("-o", "--outfolder"), type="character", default=NULL,
              help="path to output folder", metavar="character")
);

 
opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser, positional_arguments = TRUE)$options

# cat("opt: ")
# print(opt)
# cat("\n")
# cat(paste0("opt$gene_file: ", opt$gene_file, "\n"))
# cat(paste0("opt$tad_file: ", opt$tad_file, "\n"))
# cat(paste0("opt$bin_size: ", opt$bin_size, "\n"))
# cat(paste0("opt$chromo: ", opt$chromo, "\n"))
# cat(paste0("opt$outfolder: ", opt$outfolder, "\n"))

if( is.null(opt$gene_file) |is.null(opt$tad_file) |is.null(opt$bin_size) |is.null(opt$outfolder) |is.null(opt$chromo)  ){
  stop("ERROR: missing input argument \n")
}

# Rscript gene2TAD.R -f <input_file_with_genes> -t <input_file_with_TADs> -c <chromo> -o <path_to_output_folder> -b <bin_size>
# Rscript gene2TAD.R -f tmp_chr1.txt -t conservTAD_chr1_0.6.txt -c chr1 -o all_data -b 40000

chromo <- opt$chromo
binSize <- opt$bin_size
geneFile <- opt$gene_file
TADfile <- opt$tad_file
outFold <- opt$outfolder

system(paste0("mkdir -p ", outFold))

#************************************** FOR DEBUG ****************************************************************************
# chromo <- "chr1"
# binSize <- 40000
# geneFile <- "/media/electron/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt"
# TADfile <- "/media/electron/mnt/ed4/marie/TAD_call_pipeline_TopDom/consensus_TopDom_covThresh_r0.6_t80000_v0_w1_version2/chr1_conservedTADs.txt"
# outFold <- "tmp"
#*****************************************************************************************************************************

## HARD-CODED PREFIX
g2t_prefix <- paste0(outFold, "/tmp_g2t")
assigned_prefix <- paste0(outFold, "/tmp_assigned")

# geneDT <- read.delim(geneFile, header=F, stringsAsFactors = F)
# stopifnot(ncol(geneDT) == 4)
# colnames(geneDT) <- c("entrezID", "chromo", "start", "end")
# stopifnot(is.numeric(geneDT$start) & is.numeric(geneDT$end))

# x <- load(geneFile)
# geneDT <- eval(parse(text = x))
# geneDT <- na.omit(geneDT)
# geneDT <- geneDT[,c("entrezID", "Chr", "Start", "End")]
# colnames(geneDT) <- c("entrezID", "chromo", "start", "end")
# stopifnot(is.numeric(geneDT$start) & is.numeric(geneDT$end))

geneDT <- read.delim(geneFile, header=T, stringsAsFactors = F)
colnames(geneDT) <- c("entrezID", "chromo", "start", "end", "assembly", "strand")
geneDT <- geneDT[,c("entrezID", "chromo", "start", "end", "strand")]
stopifnot(is.numeric(geneDT$start) & is.numeric(geneDT$end))
geneDT <- geneDT[order(geneDT$chromo, geneDT$start),] 

#TAD_DT <- read.delim(TADfile, header=T, stringsAsFactors = F)
TAD_DT <- read.delim(TADfile, header=F, stringsAsFactors = F, col.names=c("chromo", "start", "end"))
# *** temporarily changed for Dixon2018_integrative_data & Cancer_HiC_data_TAD_DA: no header
#
#stopifnot(ncol(TAD_DT) == 5)
stopifnot(ncol(TAD_DT) == 3)
TAD_DT <- TAD_DT[,1:3]
stopifnot(all(colnames(TAD_DT) == c("chromo", "start", "end")))
stopifnot(is.numeric(TAD_DT$start) & is.numeric(TAD_DT$end))


for(chromo in unique(TAD_DT$chromo)) {

  outFile_g2t <- paste0(g2t_prefix, "_", chromo, ".txt")
  outFile_assigned <- paste0(assigned_prefix, "_", chromo, ".txt")
    
  TAD_DT_chromo <- TAD_DT[TAD_DT$chromo == chromo,]
  stopifnot(all (TAD_DT_chromo$chromo == chromo))
  stopifnot( ! any(duplicated(TAD_DT_chromo$TAD_id)))

  TAD_DT_chromo <- TAD_DT_chromo[order(TAD_DT_chromo$start),]
  
  if(TAD_DT_chromo$end[nrow(TAD_DT_chromo)] %% 10 != 0) {
    if( (TAD_DT_chromo$end[nrow(TAD_DT_chromo)] + 1) %% 10 == 0) {
	  TAD_DT_chromo$end[nrow(TAD_DT_chromo)] <- TAD_DT_chromo$end[nrow(TAD_DT_chromo)] + 1		
     } else{
		stop("even with +1 %%10 is not == 0")
	}
  }

  stopifnot( all(TAD_DT_chromo$end %% 10 == 0)  )

  # retrieve the last gene position on the current chromo
  maxEnd <- max(geneDT$end[geneDT$chromo == chromo], na.rm=T)
  
  # add a column with TAD-id
  TAD_DT_chromo$TAD_id <- paste0(chromo, "_TAD", 1:nrow(TAD_DT_chromo))
  
  # iterate over the TAD and fill with boundaries
  prev_end <- 0
  bd_count <- 1

  region_DT <- foreach(i = 1:nrow(TAD_DT_chromo), .combine = 'rbind') %do% {
    tmpDT <- data.frame(chromo = chromo, 
                        region_id = TAD_DT_chromo$TAD_id[i],
                        start = TAD_DT_chromo$start[i], 
                        end = TAD_DT_chromo$end[i]
                        )
    if(TAD_DT_chromo$start[i] > prev_end + 1 )   {
      prev_tmpDT <- data.frame(chromo = chromo, 
                               region_id = paste0(chromo, "_BOUND", bd_count),
                               start = (prev_end + 1), 
                               end = (TAD_DT_chromo$start[i] -1)
                                )
      bd_count <- bd_count + 1
      tmpDT <- rbind(prev_tmpDT, tmpDT)
    }
    prev_end <-  TAD_DT_chromo$end[i]
    tmpDT
  }
  region_DT <- as.data.frame(region_DT)
  stopifnot(ncol(region_DT) == 4)
  stopifnot( all (colnames(region_DT) == c("chromo", "region_id", "start", "end")))
  stopifnot(is.numeric(region_DT$start))
  stopifnot(is.numeric(region_DT$end))
  region_DT$region_id <- as.character(region_DT$region_id)
  region_DT$chromo <- as.character(region_DT$chromo)
  # for the last, check if the last end is > that the maxEnd (max last gene position)
  if(region_DT$end[nrow(region_DT)] < maxEnd) {
   # lastEnd <-  ceiling(249215453 / binSize) * binSize
   lastEnd <-  ceiling(maxEnd / binSize) * binSize
   cat()
    stopifnot(lastEnd >= maxEnd)  
    lastDT <- data.frame(chromo = chromo, 
                            region_id = paste0(chromo, "_BOUND", bd_count),
                             start = (region_DT$end[nrow(region_DT)] + 1), 
                             end = lastEnd
                            )
    region_DT <- rbind(region_DT, lastDT)
  }
  stopifnot( ! any(duplicated(region_DT$region_id)))
  # do some additional checks 
  # TAD id with start and end should be the same for a given TAD in both DT
  for(i in TAD_DT_chromo$TAD_id) {
    v0_line <- TAD_DT_chromo[TAD_DT_chromo$TAD_id == i, c("chromo", "TAD_id", "start", "end")]
    v1_line <- region_DT[region_DT$region_id == i,]
    stopifnot(all(v0_line == v1_line))
  }
  # the start should always be bigger than previous end
  for(i in 2:nrow(region_DT)) {
    region_DT$start[i] > region_DT$end[i-1]
  }
  # start should end with 1 and end with 0
  # and end should be a multiple of binSize, and start-1 also
  stopifnot( all(region_DT$start %% 10 == 1))  
  stopifnot( all(region_DT$end %% 10 == 0)  )
  stopifnot( all(region_DT$end %% binSize == 0)  )
  stopifnot( all( (region_DT$start-1) %% binSize == 0)  )
  
  # write the assigned regions  # paste0 to file name with chromo ??
  write.table(region_DT, file = outFile_assigned,
              sep="\t", quote=F, row.names = F, col.names = F )
  cat(paste0("... written: ", outFile_assigned, "\n"))
  
  # now assigned the genes
  # first if they have the start in a TAD region
  # then if they have the start in a TAD
  geneDT_chromo <- geneDT[geneDT$chromo == chromo, ]

  stopifnot(is.numeric(geneDT_chromo$start[1]))
  stopifnot(is.numeric(geneDT_chromo$end[1]))
  
  gene2tad_DT <- foreach(i = 1:nrow(geneDT_chromo), .combine = 'rbind') %do% {
    gene_start <- geneDT_chromo$start[i]
    gene_end <- geneDT_chromo$end[i]
    gene_strand <- as.character(geneDT_chromo$strand[i])
    # in which region we find the start
    region_start <- region_DT$region_id[region_DT$start <= gene_start & region_DT$end >= gene_start ]
    stopifnot(length(region_start) > 0)
    # in which region we find the end
    region_end <- region_DT$region_id[region_DT$start <= gene_end & region_DT$end >= gene_end ]
    stopifnot(length(region_end) > 0)
    
    #### UPDATE 22.09 => IF THE GENE IS ON THE NEGATIVE STRAND -> THE START IS IN FACT THE END
    if(gene_strand == "-") {
      tmpvar <- region_start
      region_start <- region_end
      region_end <- tmpvar
    }
    
    # if the start is a TAD -> TAD of the start
    # if the start is not in a TAD but the end is -> TAD of the end
    # else region of the start
    gene_region <- ifelse( regexpr("_TAD", region_start) > 0, region_start,
                           ifelse(regexpr("_TAD", region_end) > 0, region_end, region_start))
                           
    data.frame(entrezID = geneDT_chromo$entrezID[i],
               chromo = geneDT_chromo$chromo[i],
               start = geneDT_chromo$start[i],
               end = geneDT_chromo$end[i],
                        region = gene_region)
  }
  stopifnot(ncol(gene2tad_DT) == 5)
  stopifnot(nrow(gene2tad_DT) == nrow(geneDT_chromo))
  stopifnot(all(colnames(gene2tad_DT) == c("entrezID", "chromo", "start", "end", "region")))
  gene2tad_DT$entrezID <- as.character(gene2tad_DT$entrezID)
  gene2tad_DT$region <- as.character(gene2tad_DT$region)
  stopifnot(all(regexpr("_TAD", gene2tad_DT$region) > 0 | regexpr("_BOUND", gene2tad_DT$region) > 0))
  stopifnot(all(gene2tad_DT$chromo == gsub("(.+)_.+", "\\1", gene2tad_DT$region)))
  
  # write the assigned regions  # paste0 to file name with chromo ??
  write.table(gene2tad_DT, file = outFile_g2t,
              sep="\t", quote=F, row.names = F, col.names = F )
  cat(paste0("... written: ", outFile_g2t, "\n"))
  
}



