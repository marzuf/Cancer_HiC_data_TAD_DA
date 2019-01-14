#!/usr/bin/bash

# ./aws2_all_chromo_hicData_rebinning.sh K562
# (see what has been run in all_cmds_pip.txt)

start_time=$(date -R)    

if [[ $# != 1 ]]; then
    echo "invalid # of arguments"
    exit 1
fi
##################### FUNCTION DEFINITION
checkFile() {
if [[ ! -f  $2 ]]; then
    echo "... $1 ($2) does not exist !"
    exit 1
fi
}
##########################################
cl="$1"

rexec="Rscript"

script_rebinning="aws2_prepFromHiC_aggregCounts.R"
checkFile script_rebinning $script_rebinning

script_check="cmp_rebinning_custom_HiTC.R"
checkFile script_check $script_check

oldBinSizeKb="10"
newBinSizeKb="40"

#all_chromo=( "chr"{1..22} "chrX" )
#all_chromo=( "chr1" "chr9" )
all_chromo=( "chrX" )

echo "*** START ***"
echo "... > Cell line: $cl"
echo "... > Chromosome(s): ${all_chromo[*]}"

oldFold="${cl}/RAW_${oldBinSizeKb}kb"

newFold="${cl}/RAW_${newBinSizeKb}kb"
mkdir -p $newFold

for chromo in "${all_chromo[@]}"; do
	
#  Rscript prepFromHiC_aggregCounts.R -f MCF-7/RAW_10kb/MCF-7_chr16_RAW_10kb.hic.counts -F MCF-7/RAW_40kb/MCF-7_chr16_RAW_40kb.hic.counts -b 10000 -B 40000 > MCF-7/RAW_40kb/MCF-7_chr16_RAW_40kb.hic.counts_aggregCounts.log
#  Rscript prepFromHiC_aggregCounts.R -f oldFile -F newFile -b oldBinSizeBp -B newBinSizeBp > logfile

	oldCountFile="$oldFold/${cl}_${chromo}_RAW_${oldBinSizeKb}kb.hic.counts"
	newCountFile="$newFold/${cl}_${chromo}_RAW_${newBinSizeKb}kb.hic.counts"
	
	checkFile oldCountFile $oldCountFile

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
