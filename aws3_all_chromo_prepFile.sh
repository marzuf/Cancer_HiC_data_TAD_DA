#!/usr/bin/bash

# ./aws3_all_chromo_prepFile.sh K562
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

script_pre="aws3_prepFromHiC_prepFile.R"
checkFile script_pre $script_pre

binSizeKb="40"

all_chromo=( "chr"{1..22} "chrX" )
#all_chromo=( "chr21" )
#all_chromo=( "chr1" "chr9" )

inFold="${cl}/RAW_${binSizeKb}kb"

outFold="${cl}/PRE_${binSizeKb}kb"
mkdir -p $outFold

for chromo in "${all_chromo[@]}"; do
	#Rscript prepFromHiC_prepFile.R -f MCF-7/RAW_40kb/MCF-7_chr21_RAW_40kb.hic.counts -F MCF-7/PRE_40kb/MCF-7_chr21_RAW_40kb.pre -c chr21 -b 40000
	#Rscript prepFromHiC_prepFile.R -f inFile_aggregCounts -F outFile_pre -c chromo -b binSize

	inFile="$inFold/${cl}_${chromo}_RAW_${binSizeKb}kb.hic.counts"
	checkFile inFile $inFile
	
	outFile="$outFold/${cl}_${chromo}_RAW_${binSizeKb}kb.pre"

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
