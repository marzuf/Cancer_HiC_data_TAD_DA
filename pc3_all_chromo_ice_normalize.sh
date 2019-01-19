#!/usr/bin/bash

# ./pc3_all_chromo_ice_normalize.sh GSE73782_PC3
# (see what has been run in all_cmds_pip.txt)

if [[ $# != 1 ]]; then
    echo "invalid # of arguments"
    exit 1
fi

cl="$1"

start_time=$(date -R)    

ice_norm_rscript="rao_ice_normalize.R"
resolKb="40"

inFold="$cl/RAW_${resolKb}kb"
outFold="${cl}_${resolKb}kb/TopDom_MAT"
mkdir -p $outFold



all_chromo=( "chr"{1..22} "chrX" )
#all_chromo=( "1" "6" "15" "21" )
#all_chromo=( "chr20" )

echo "*** START ***"
echo "... > Cell line: $cl"
echo "... > Chromosome(s): ${all_chromo[*]}"


for chromo in "${all_chromo[@]}"; do

	cmd="Rscript $ice_norm_rscript $inFold/${cl}_${chromo}_RAW_${resolKb}kb.hic.counts $outFold/${cl}_${chromo}_ICE_${resolKb}kb.hic.TopDom.matrix $chromo ${resolKb}000"

	echo $cmd
	$cmd


done

########## END ####################################################################################################################################
###################################################################################################################################################
echo "*** DONE"
end_time=$(date -R)    
echo $start_time
echo $end_time
exit 0
