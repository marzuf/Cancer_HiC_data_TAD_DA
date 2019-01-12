#!/usr/bin/bash

# ./all_chromo_prepFile.sh

start_time=$(date -R)    

rexec="Rscript"

script_pre="prepFromHiC_prepFile.R"

binSizeKb="40"

all_chromo=( "chr"{1..22} "chrX" )
#all_chromo=( "chr21" )

inFold="MCF-7/RAW_${binSizeKb}kb"

outFold="MCF-7/PRE_${binSizeKb}kb"
mkdir -p $outFold

for chromo in "${all_chromo[@]}"; do
	#Rscript prepFromHiC_prepFile.R -f MCF-7/RAW_40kb/MC-7_chr21_RAW_40kb.hic.counts -F MCF-7/PRE_40kb/MC-7_chr21_RAW_40kb.pre -c chr21 -b 40000
	#Rscript prepFromHiC_prepFile.R -f inFile_aggregCounts -F outFile_pre -c chromo -b binSize

	inFile="$inFold/MC-7_${chromo}_RAW_${binSizeKb}kb.hic.counts"
	if [[ ! -f $inFile ]]; then
		echo "!!! inFile: $inFile does not exists !!!"
		exit 1
	fi
	
	outFile="$outFold/MC-7_${chromo}_RAW_${binSizeKb}kb.pre"

	cmd_pre="$rexec $script_pre -f $inFile -F $outFile -c $chromo -b ${binSizeKb}000" 
	echo $cmd_pre
	$cmd_pre

done

########## END ####################################################################################################################################
###################################################################################################################################################
echo "*** DONE"
end_time=$(date -R)    
echo $start_time
echo $end_time
exit 0
