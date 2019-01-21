#!/usr/bin/bash

# ./1_all_chromo_hic_to_TADs.sh MCF-7  # start from step2 if download from juicer AWS - 14.01.2018
# ./1_all_chromo_hic_to_TADs.sh ENCSR549MGQ_T47D # - 14.01.2018
# ./1_all_chromo_hic_to_TADs.sh K562 # start from step2 if download from juicer AWS - 14.01.2018
# ./1_all_chromo_hic_to_TADs.sh ENCSR444WCZ_A549 # - 14.01.2018
# (see what has been run in all_cmds_pip.txt)

# if not downloaded from Juicer, start with hic data located in 
# "/mnt/ndata/Yuanlong/2.Results/1.Juicer/mega_<cell_line>/mega/aligned/inter_30.hic"

start_time=$(date -R)    
# set -e                #  in case there is no chr23, otherwise not run chrX

if [[ $# != 1 ]]; then
    echo "invalid # of arguments"
    exit 1
fi


clName="$1"
#chromo="$2"

all_chromo=( "chr"{1..22} "chrX" )
#all_chromo=( "chr15" "chr16" "chr17" )
#all_chromo=( "chr1" "chr9" )
#all_chromo=( "chrX" )

#all_chromo=( "chr21" )


echo "*** START ***"
echo "... > Cell line: $clName"
echo "... > Chromosome(s): ${all_chromo[*]}"

step1=1		# extract chromo matrix from hic file # <cell_line>/NORM_MAT
step2=1		# convert Rao format to TopDom format # <cell_line>/TopDom_MAT
step3=1		# run TopDom # <cell_line>/TopDom_FILES
step4=1     # convert TopDom BED format to FINAL DOMAINS BED format # <cell_line>/FINAL_DOMAINS

# INPUT: .hic files from merged replicates

##################################################**** SOME FUNCTIONS
runCMD() {
  echo "> $1"
  eval $1
}

checkFile() {
if [[ ! -f  $2 ]]; then
    echo "... $1 ($2) does not exist !"
    exit 1
fi
}
#############################################************************


#####**** HARD-CODED SETTINGS
# for all steps
binSizeKb="40"
Rexec="Rscript"

	if [[ "$clName" ==  "GSE73782_PC3_ICE" ]]; then
		norm="ICE"
	else
		norm="KR"
	fi


#norm="KR"

for chromo in "${all_chromo[@]}"; do

	# step1:
	step1_outFolder="NORM_MAT"
	juicerExec="java -Xms512m -Xmx2048m -jar /mnt/ed2/shared/TADcompare/Software/juicer/juicer_tools.jar"
	step1_outFile="${clName}_${binSizeKb}kb/$step1_outFolder/${clName}_${chromo}_${norm}_${binSizeKb}kb.hic.matrix"
	hicFolder="/mnt/ndata/Yuanlong/2.Results/1.Juicer"
	
	
	if [[ "$clName" ==  "GSE105318_DLD1" ]]; then
		hicFile="$hicFolder/mega_$clName/aligned/inter_30.hic"
	else
		hicFile="$hicFolder/mega_$clName/mega/aligned/inter_30.hic"
	fi
	

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

	checkFile hicFile $hicFile

	runCMD "mkdir -p `dirname $step1_outFile`"
	runCMD "$juicerExec dump observed $norm $hicFile $chromo $chromo BP ${binSizeKb}000 $step1_outFile"

	checkFile step1_outFile $step1_outFile


	fi

	# write: PC3_rep_12_100kb/NORM_MAT/PC3_rep_12_chr1_KR_100kb.hic.matrix

	if [[ $step2 -eq 1 ]]; then
	
	######################################################################################
	# STEP 2: convert Rao normalized matrix to TopDom format
	######################################################################################
	# Rscript convert_Rao_TopDom.R inFile outFile chromo binSizeBp
	# save file in $clName/TopDom_MAT

	checkFile step1_outFile $step1_outFile

	runCMD "$Rexec $step2_rao2td_script $step1_outFile $step2_outFile $chromo ${binSizeKb}000"
	
	checkFile step2_outFile $step2_outFile

	fi



	# written: PC3_rep_12_100kb/TopDom_MAT/PC3_rep_12_chr1_KR_100kb.hic.TopDom.matrix

	if [[ $step3 -eq 1 ]]; then
	######################################################################################
	# STEP 3: run TopDom
	######################################################################################
	# Rscript convert_Rao_TopDom.R inFile outFile chromo binSizeBp
	# save file in $clName/TopDom_FILES

	checkFile step2_outFile $step2_outFile

	runCMD "$Rexec $step3_topdom_script -i $step2_outFile -o $step3_outPrefix -w $window_size"

	checkFile step3_outFile $step3_outFile


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

	checkFile step3_outFile $step3_outFile

	runCMD "mkdir -p `dirname $step4_outFile`"
	runCMD "grep domain $step3_outFile | cut -f1-3 | awk '{print \$1 \"\t\" (\$2 + 1) \"\t\" \$3}' > $step4_outFile"


	checkFile step4_outFile $step4_outFile

	runCMD "head -2 $step4_outFile"

	fi

done # end iterating over chromosomes

###################################################################################################################################################
########## END ####################################################################################################################################
###################################################################################################################################################
echo "*** DONE"
end_time=$(date -R)    
echo $start_time
echo $end_time
exit 0










































