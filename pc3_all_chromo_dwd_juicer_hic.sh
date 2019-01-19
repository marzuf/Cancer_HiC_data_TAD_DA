#!/usr/bin/bash

# ./pc3_all_chromo_dwd_juicer_hic.sh GSE73782_PC3
# (see what has been run in all_cmds_pip.txt)

if [[ $# != 1 ]]; then
    echo "invalid # of arguments"
    exit 1
fi

cl="$1"

if [[ $cl == "GSE73782_PC3" ]]; then
	url="/mnt/ndata/Yuanlong/2.Results/1.Juicer/mega_${cl}/mega/aligned/inter_30.hic"
fi

start_time=$(date -R)    

javaExec="java -Xms512m -Xmx2048m -jar /mnt/ed2/shared/TADcompare/Software/juicer/juicer_tools.jar"

resolKb="40"

all_chromo=( {1..22} "X" )
#all_chromo=( "1" "6" "15" "21" )
#all_chromo=( "21" )

echo "*** START ***"
echo "... > Cell line: $cl"
echo "... > Chromosome(s): ${all_chromo[*]}"


outFold="$cl/RAW_${resolKb}kb"
mkdir -p $outFold

for chromo in "${all_chromo[@]}"; do

	outFile="$outFold/${cl}_chr${chromo}_RAW_${resolKb}kb.hic.counts"

	cmd="$javaExec dump observed None $url $chromo $chromo BP ${resolKb}000 $outFile"

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
