#!/usr/bin/bash

# ./3_assign_genes.sh

start_time=$(date -R)    
set -e

#if [[ $# != 2 ]]; then
if [[ $# != 1 ]]; then

    echo "invalid # of arguments"
    exit 1
fi


clName="$1"
#chromo="$2"

echo "*** START ***"
echo "... > Cell line: $clName"
#echo "... > Chromosome: $chromo"

step1=1		# assign genes to TADs


# input file:
# MCF-7_40kb/FINAL_DOMAINS/MCF-7_chr15_KR_40kb_final_domains.txt
#  <cl1cl2>_40kb/FINAL_DOMAINS/<cl1cl2>_chr15_KR_40kb_final_domains.txt  # for the consensus: pass cl1namecl2name as clname
 
# .../${clName}_${chromo}_${norm}_${binSize}kb_final_domains.txt
  
# output file:

#####**** SOME FUNCTIONS

runCMD() {
  echo "$1"
  eval $1
}
export -f runCMD

checkFile() {
if [[ ! -f  $2 ]]; then
    echo "... $1 ($2) does not exist !"
    exit 1
fi
}

#####************************


#####**** HARD-CODED SETTINGS

# for all steps
norm="KR"
binSizeKb="40"
Rexec="Rscript"
mainFold="${clName}_${binSizeKb}kb"

maxJobs=40
maxLoad=60

#all_chromo=( "chr"{1..22} "chrX" )  # -> need to do all chromo to have 1 single file at the end
all_chromo=( "chr15" )

# step1:
g2t_script="gene2TAD_consensus_version2.R"
infold_genes="/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt"
# HARD CODED IN g2t assignment R SCRIPT (used at the end  of STEP 1 for concatenating files)
g2t_prefix="tmp_g2t"
reg_prefix="tmp_assigned"

step1_inputFolder="$mainFold/FINAL_DOMAINS"

step1_outFolder_tmp="$mainFold/genes2tad/tmp"
step1_outFolder_final="$mainFold/genes2tad"

step1_outFile_genes="$step1_outFolder_final/all_genes_positions.txt"
step1_outFile_tads="$step1_outFolder_final/all_assigned_regions.txt"

#####************************

if [[ $step1 -eq 1 ]]; then

runCMD "mkdir -p $step1_outFolder_final"
runCMD "mkdir -p $step1_outFolder_tmp"

# clean tmp folder (cat command at the end !)
if [[ -d $step1_outFolder_tmp ]]; then
echo "delete folder"
	runCMD "rm -r $step1_outFolder_tmp"
fi

parallel -i -j $maxJobs -l $maxLoad sh -c "runCMD 'Rscript $g2t_script -f $infold_genes -t $step1_inputFolder/${clName}_{}_${norm}_${binSizeKb}kb_final_domains.txt -c {} -o $step1_outFolder_tmp -b ${binSizeKb}000'" -- ${all_chromo[@]}

runCMD "cat $step1_outFolder_tmp/$g2t_prefix* > $step1_outFile_genes"
runCMD "cat $step1_outFolder_tmp/$reg_prefix* > $step1_outFile_tads"

checkFile step1_outFile_genes $step1_outFile_genes
checkFile step1_outFile_tads $step1_outFile_tads


fi

###################################################################################################################################################
########## END ####################################################################################################################################
###################################################################################################################################################

echo "*** DONE"
end_time=$(date -R)    
echo $start_time
echo $end_time


exit 0