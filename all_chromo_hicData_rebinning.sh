#!/usr/bin/bash

# ./all_chromo_hicData_rebinning.sh

start_time=$(date -R)    

rexec="Rscript"

script_rebinning="prepFromHiC_aggregCounts.R"

script_check="cmp_rebinning_custom_HiTC.R"

oldBinSizeKb="10"
newBinSizeKb="40"

all_chromo=( "chr"{1..22} "chrX" )
#all_chromo=( "chr21" )

oldFold="MCF-7/RAW_${oldBinSizeKb}kb"

newFold="MCF-7/RAW_${newBinSizeKb}kb"

mkdir -p $newFold

for chromo in "${all_chromo[@]}"; do

	
#  Rscript prepFromHiC_aggregCounts.R -f MCF-7/RAW_10kb/MCF-7_chr16_RAW_10kb.hic.counts -F MCF-7/RAW_40kb/MCF-7_chr16_RAW_40kb.hic.counts -b 10000 -B 40000 > MCF-7/RAW_40kb/MCF-7_chr16_RAW_40kb.hic.counts_aggregCounts.log
#  Rscript prepFromHiC_aggregCounts.R -f oldFile -F newFile -b oldBinSizeBp -B newBinSizeBp > logfile


	oldCountFile="$oldFold/MCF-7_${chromo}_RAW_${oldBinSizeKb}kb.hic.counts"
	newCountFile="$newFold/MCF-7_${chromo}_RAW_${newBinSizeKb}kb.hic.counts"
	
	if [[ ! -f $oldCountFile ]]; then
		echo "!!! oldCountFile: $oldCountFile does not exists !!!"
		exit 1
	fi



	cmd_rebin="$rexec $script_rebinning -f $oldCountFile -F $newCountFile -b ${oldBinSizeKb}000 -B ${newBinSizeKb}000" # will not work with flags> ${newCountFile}_aggregCounts.log"
	echo $cmd_rebin
	$cmd_rebin
	
	
# Rscript cmp_rebinning_custom_HiTC.R MCF-7/RAW_10kb/MC-7_chr16_RAW_10kb.hic.counts MCF-7/RAW_40kb/MCF-7_chr16_RAW_40kb.hic.counts 10000 40000 chr16 > MCF-7/RAW_40kb/MCF-7_chr16_RAW_40kb.hic.counts_checkRebinning.log
# Rscript cmp_rebinning_custom_HiTC.R inFile outFile oldBinBp newBinBp chromo > logfile

	
	
	cmd_check="$rexec $script_check $oldCountFile $newCountFile ${oldBinSizeKb}000 ${newBinSizeKb}000 $chromo" #> ${newCountFile}_checkRebinning.log"
	echo $cmd_check
	$cmd_check


done

########## END ####################################################################################################################################
###################################################################################################################################################
echo "*** DONE"
end_time=$(date -R)    
echo $start_time
echo $end_time
exit 0
