#!/usr/bin/bash

# ./prep_assign_genes_singleDS.sh leukemia/K562/GSE63525/TopDom GSE63525_K562

# ./prep_assign_genes_singleDS.sh colon/DLD1/GSE105318/ENCFF439QFU/TopDom GSE105318_ENCFF439QFU_DLD1

# ./prep_assign_genes_singleDS.sh pancreas/Panc1/GSE105566/ENCFF358MNA/TopDom GSE105566_ENCFF358MNA_Panc1

#set -e


if [[ $# != 2 ]]; then
    echo "invalid # of arguments"
    exit 1
fi

start_time=$(date -R)    

Rexec="Rscript"

## PARALLELIZATION
maxJobs=10
maxLoad=70

# general settings
chromo=( "chr"{1..22} "chrX" ) 
#chromo=( "chr21" ) 
#chromo=( "chr1" )
ncpu="5"

inputFolder="$1"
outDatasetName="$2"

# take the TopDom TAD files -> put them in correct format in FIND_CONSENSUS_TADS => so that then I can run assign_genes.sh for single dataset TADs (not only consensus TADs)
#head leukemia/K562/GSE63525/TopDom/GSE63525_K562_40kb_ICE_chr16_final_domains.txt
#chr16   80001   440000
#chr16   440001  760000

# e.g. the format should look like
#head FIND_CONSENSUS_TADS/GSE105194_ENCFF027IEO_astroCerebellum_vs_GSE105957_ENCFF715HDW_astroSpinal/chr1_conservedTADs.txt 
#chr1    1       680000
#chr1    680001  1520000

# e.g. transfer
# leukemia/K562/GSE63525/TopDom/GSE63525_K562_40kb_ICE_chr16_final_domains.txt
# to
# FIND_CONSENSUS_TADS/GSE63525_K562/chr16_conservedTADs.txt 

outFolder="/mnt/etemp/marie/Dixon2018_integrative_data/FIND_CONSENSUS_TADS/$outDatasetName"

echo "> mkdir -p $outFolder"
mkdir -p $outFolder

parallel -i -j $maxJobs -l $maxLoad sh -c "echo cp $inputFolder/*_{}_final_domains.txt $outFolder/{}_conservedTADs.txt" -- ${chromo[@]}
parallel -i -j $maxJobs -l $maxLoad sh -c "cp $inputFolder/*_{}_final_domains.txt $outFolder/{}_conservedTADs.txt" -- ${chromo[@]}

echo "> ls -1 $outFolder"
ls -1 $outFolder


###################################################################################################################################################
########## END ####################################################################################################################################
###################################################################################################################################################

echo "*** DONE"
end_time=$(date -R)    
echo $start_time
echo $end_time


