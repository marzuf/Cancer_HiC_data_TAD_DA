#!/usr/bin/bash

# ./all_chromo_pre2hic2matrix.sh

start_time=$(date -R)    

javaExec="java -Xms512m -Xmx2048m -jar /mnt/ed2/shared/TADcompare/Software/juicer/juicer_tools.jar"

resolKb="40"

all_chromo=( {1..22} "X" )
#all_chromo="21"

inFold="MCF-7/PRE_${resolKb}kb"

outFold="MCF-7/RAW_HIC_INTRA_${resolKb}kb"
mkdir -p $outFold

step1_outFolder="MCF-7_${resolKb}kb/NORM_MAT"
mkdir -p $step1_outFolder

norm="KR"


for chromo in "${all_chromo[@]}"; do

	inFile="$inFold/MCF-7_chr${chromo}_RAW_${resolKb}kb.pre"
	if [[ ! -f  $inFile ]]; then
   	 echo "... inFile ($inFile) does not exist !"
   	 exit 1
	fi
	
	sizeFile="$inFold/chr${chromo}.size"
	if [[ ! -f  $sizeFile ]]; then
   	 echo "... sizeFile ($sizeFile) does not exist !"
   	 exit 1
	fi

	outFile="$outFold/MCF-7_chr${chromo}_intra_${resolKb}kb.hic"

	cmd="$javaExec pre -r ${resolKb}000 -c chr$chromo $inFile $outFile $sizeFile"
	echo $cmd
	$cmd
	
	step1_outFile="$step1_outFolder/MCF-7_chr${chromo}_${norm}_${resolKb}kb.hic.matrix"
	cmdNorm="$javaExec dump observed $norm $outFile $chromo $chromo BP ${resolKb}000 $step1_outFile"
	echo $cmdNorm
	$cmdNorm

done

########## END ####################################################################################################################################
###################################################################################################################################################
echo "*** DONE"
end_time=$(date -R)    
echo $start_time
echo $end_time
exit 0

# java -Xms512m -Xmx2048m -jar /mnt/ed2/shared/TADcompare/Software/juicer/juicer_tools.jar pre -r 40000 -d -c 21 MCF-7/PRE_40kb/MCF-7_chr21_RAW_40kb.pre MCF-7/HIC_INTRA_30_40kb/MCF-7_chr21_intra_40kb.hic hg19

# java -Xms512m -Xmx2048m -jar /mnt/ed2/shared/TADcompare/Software/juicer/juicer_tools.jar pre -r 50000 -d -c 6 /mnt/etemp/marie/NAS_FOLDERS/TADcompare_benchmark/pipeline_chr6/12_07_Cmap_50k_rep1_ICE/input_caller/arrowhead_Cmap_mapq30_diffrag_intrachr6_bw50kb_subsample_rep3_ICE_pre.txt foo_hic_chr6.hic hg19

#    !!! -d => only intrachromo !!!    
#    pre -q mapqScore -r resolution -d -c chromo infile outfile hg19
#pre -q 30 -r 40000 -d -c 21 infile outfile hg19
#MCF-7/PRE_40kb/MCF-7_chr21_RAW_40kb.pre


