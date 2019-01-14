#!/usr/bin/bash

# ./aws4_all_chromo_pre2hic2matrix.sh K562
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

javaExec="java -Xms512m -Xmx2048m -jar /mnt/ed2/shared/TADcompare/Software/juicer/juicer_tools.jar"

resolKb="40"

all_chromo=( {1..22} "X" )
# all_chromo=( "1" "9" )

echo "*** START ***"
echo "... > Cell line: $cl"
echo "... > Chromosome(s): ${all_chromo[*]}"

inFold="${cl}/PRE_${resolKb}kb"

outFold="${cl}/RAW_HIC_INTRA_${resolKb}kb"
mkdir -p $outFold

step1_outFolder="${cl}_${resolKb}kb/NORM_MAT"
mkdir -p $step1_outFolder

norm="KR"

for chromo in "${all_chromo[@]}"; do

	inFile="$inFold/${cl}_chr${chromo}_RAW_${resolKb}kb.pre"
	if [[ ! -f  $inFile ]]; then
   	 echo "... inFile ($inFile) does not exist !"
   	 exit 1
	fi
	
	sizeFile="$inFold/chr${chromo}.size"
	if [[ ! -f  $sizeFile ]]; then
   	 echo "... sizeFile ($sizeFile) does not exist !"
   	 exit 1
	fi

	outFile="$outFold/${cl}_chr${chromo}_intra_${resolKb}kb.hic"

	cmd="$javaExec pre -r ${resolKb}000 -c chr$chromo $inFile $outFile $sizeFile"
	echo $cmd
	$cmd
	
	step1_outFile="$step1_outFolder/${cl}_chr${chromo}_${norm}_${resolKb}kb.hic.matrix"
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

