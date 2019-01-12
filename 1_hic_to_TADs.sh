#!/usr/bin/bash

# ./1_hic_to_TADs.sh PC3_rep12 chr1
# ./1_hic_to_TADs.sh MC-7 chr15


start_time=$(date -R)    
set -e

if [[ $# != 2 ]]; then
    echo "invalid # of arguments"
    exit 1
fi


clName="$1"
chromo="$2"

echo "*** START ***"
echo "... > Cell line: $clName"
echo "... > Chromosome: $chromo"

step1=0		# extract chromo matrix from hic file
step2=0		# convert Rao format to TopDom format
step3=1		# run TopDom
step4=0     # convert TopDom BED format to FINAL DOMAINS BED FORMAT

# INPUT: .hic files from merged replicates

#####**** SOME FUNCTIONS

function runCMD {
  echo "> $cmd"	
  $cmd  
}

#####************************


#####**** HARD-CODED SETTINGS

# for all steps
norm="KR"
binSizeKb="40"
Rexec="Rscript"

# step1:
step1_outFolder="NORM_MAT"
juicerExec="java -Xms512m -Xmx2048m -jar /mnt/ed2/shared/TADcompare/Software/juicer/juicer_tools.jar"
step1_outFile="${clName}_${binSizeKb}kb/$step1_outFolder/${clName}_${chromo}_${norm}_${binSizeKb}kb.hic.matrix"
hicFolder="/mnt/ndata/Yuanlong/2.Results/1.Juicer"
hicFile="$hicFolder/mega_$clName/mega/aligned/inter_30.hic"

# step2:
step2_outFolder="TopDom_MAT"
step2_rao2td_script="convert_Rao_TopDom.R"
step2_outFile="${clName}_${binSizeKb}kb/$step2_outFolder/${clName}_${chromo}_${norm}_${binSizeKb}kb.hic.TopDom.matrix"

# step3:
step3_outFolder="TopDom_FILES"
step3_topdom_script="/mnt/ed2/shared/TADcompare/Software/TopDom/run_TopDom.R"
window_size=5
step3_outPrefix="${clName}_${binSizeKb}kb/$step3_outFolder/${clName}_${chromo}_${norm}_${binSizeKb}kb"

# step4:
step3_outFile="${step3_outPrefix}.bed"
step4_outFolder="FINAL_DOMAINS"
step4_outFile="${clName}_${binSizeKb}kb/$step4_outFolder/${clName}_${chromo}_${norm}_${binSizeKb}kb_final_domains.txt"

#####************************


#####**** RETRIEVED FROM COMMAND LINE
#clName="PC3_rep12"
#chromo="chr1"

#####************************

if [[ $step1 -eq 1 ]]; then

######################################################################################
# STEP 1: extract intrachromo 40 kb matrices from .hic files
######################################################################################
# java -Xms512m -Xmx2048m -jar juicer_tools.jar dump observed KR inFile chromo chromo BP binSizeKb outFile
# save file in $clName/NORM_MAT
if [[ ! -f  $hicFile ]]; then
    echo "... hicFile ($hicFile) does not exist !"
    exit 1
fi
cmd="mkdir -p `dirname $step1_outFile`"
runCMD $cmd

cmd="$juicerExec dump observed $norm $hicFile $chromo $chromo BP ${binSizeKb}000 $step1_outFile"
runCMD $cmd
fi

# write: PC3_rep_12_100kb/NORM_MAT/PC3_rep_12_chr1_KR_100kb.hic.matrix

if [[ $step2 -eq 1 ]]; then
######################################################################################
# STEP 2: convert Rao normalized matrix to TopDom format
######################################################################################
# Rscript convert_Rao_TopDom.R inFile outFile chromo binSizeBp
# save file in $clName/TopDom_MAT
if [[ ! -f  $step1_outFile ]]; then
    echo "... step1_outFile ($step1_outFile) does not exist !"
    exit 1
fi
cmd="$Rexec $step2_rao2td_script $step1_outFile $step2_outFile $chromo ${binSizeKb}000"
runCMD $cmd
fi

# written: PC3_rep_12_100kb/TopDom_MAT/PC3_rep_12_chr1_KR_100kb.hic.TopDom.matrix

if [[ $step3 -eq 1 ]]; then
######################################################################################
# STEP 3: run TopDom
######################################################################################
# Rscript convert_Rao_TopDom.R inFile outFile chromo binSizeBp
# save file in $clName/TopDom_FILES
if [[ ! -f  $step2_outFile ]]; then
    echo "... step2_outFile ($step2_outFile) does not exist !"
    exit 1
fi
cmd="$Rexec $step3_topdom_script -i $step2_outFile -o $step3_outPrefix -w $window_size"
runCMD $cmd
fi

# Rscript /mnt/ed2/shared/TADcompare/Software/TopDom/run_TopDom.R -i GM12878_50kb/TopDom_MAT/GM12878_chr1_KR_50kb.hic.TopDom.matrix -o GM12878_50kb/TopDom_FILES/GM12878_chr1_KR_50kb -w 5
# GM12878_50kb/TopDom_FILES/GM12878_chr1_KR_50kb.bed
# step3_outPrefix: GM12878_50kb/TopDom_FILES/GM12878_chr1_KR_50kb

if [[ $step4 -eq 1 ]]; then
######################################################################################
# STEP 4: convert TopDom BED format to "final_domains.txt" format
######################################################################################
# grep domain inFile | cut -f1-3 > outFile
# save file in $clName/FINAL_DOMAINS
if [[ ! -f  $step3_outFile ]]; then
    echo "... step3_outFile ($step3_outFile) does not exist !"
    exit 1
fi
cmd="grep domain $step3_outFile | cut -f1-3 > $step4_outFile"
runCMD $cmd
fi

###################################################################################################################################################
########## END ####################################################################################################################################
###################################################################################################################################################
echo "*** DONE"
end_time=$(date -R)    
echo $start_time
echo $end_time
exit 0










































