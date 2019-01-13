#!/usr/bin/bash

# ./all_chromo_dwd_juicer_hic_MCF7.sh

start_time=$(date -R)    

javaExec="java -Xms512m -Xmx2048m -jar /mnt/ed2/shared/TADcompare/Software/juicer/juicer_tools.jar"

resolKb="10"

#all_chromo=( {1..22} "X" )
all_chromo=( "1" "6" "15" "21" )


outFold="MCF-7/RAW_${resolKb}kb"
mkdir -p $outFold

for chromo in "${all_chromo[@]}"; do

	outFile="$outFold/MCF-7_chr${chromo}_RAW_${resolKb}kb.hic.counts"

	cmd="$javaExec dump observed None https://hicfiles.s3.amazonaws.com/external/barutcu/MCF-7.hic $chromo $chromo BP ${resolKb}000 $outFile"

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
