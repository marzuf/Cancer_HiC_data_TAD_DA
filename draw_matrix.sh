#!/usr/bin/bash

# ./draw_matrix.sh

#Rscript draw_matrix.R \
#-m /mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_aggregMaps/AGGREG_MAPS_RESCALE_TOTCOUNT_NORMmat/chr21/chr21_aggregMap.txt \
#-o DRAW_METATAD/RESCALE_TOTCOUNT_NORMmat/chr21_draw_test_meanInterCount_hierarchLev25_35280001_40560000.png \
#-k 0 \
#-b 40000 \
#-s 35280001 \
#-e 40560000 \
#-c chr21 \
#-d 25

# !!! warning hard-coded settings in draw_matrix_metaTAD.R !!!
#- matrixHeader 
#- tadHeader
#- metaTADheader 
#- featureHeader 

maxJobs=100
maxLoad=100

script_name="draw_matrix.R"

binSizeKb=40
colToSkip=3

all_chromos=( "chr1" "chr9" "chr21" )
#all_chromos=( "chr5" )
#all_chromos=( "chr"{1..22} "chrX" ) 

all_dataset_id=(
"ENCSR549MGQ_T47D"
"K562"
"ENCSR444WCZ_A549"
"NCI-H460"
"ENCSR346DCU_LNCaP"
"Panc1_rep12"
"ENCSR079VIJ_G401"
"ENCSR401TBQ_Caki2"
"ENCSR312KHQ_SK-MEL-5"
"ENCSR862OGI_RPMI-7951"
"GSM2334834_U266_HindIII"
"GSM2334832_RPMI-8226_HindIII"
"ENCSR834DXR_SK-N-MC"
"ENCSR105KFX_SK-N-DZ"
)


for dataID in "${all_dataset_id[@]}"; do


	matrixFolder="${dataID}_${binSizeKb}kb/TopDom_MAT"
	matrixPrefix="${dataID}_"
	matrixSuffix="_KR_${binSizeKb}kb.hic.TopDom.matrix"
	outFolder="DRAW_ENCODE_MATRIX/${dataID}_${binSizeKb}kb"



    parallel -i -j $maxJobs -l $maxLoad sh -c "echo Rscript $script_name  -m $matrixFolder/$matrixPrefix{}$matrixSuffix -o $outFolder/{}_$matrixPrefix{}${matrixSuffix}.png -k $colToSkip -b ${binSizeKb}000 -c {}" -- ${all_chromos[@]}
    parallel -i -j $maxJobs -l $maxLoad sh -c "Rscript $script_name  -m $matrixFolder/$matrixPrefix{}$matrixSuffix -o $outFolder/{}_$matrixPrefix{}${matrixSuffix}.png -k $colToSkip -b ${binSizeKb}000 -c {}" -- ${all_chromos[@]}


    start_positions=( "35280001" )
    end_positions=( "40560000" )

    for i in "${!start_positions[@]}"; do
        start_pos="${start_positions[i]}"
        end_pos="${end_positions[i]}"

        parallel -i -j $maxJobs -l $maxLoad sh -c "echo Rscript $script_name  -m $matrixFolder/$matrixPrefix{}$matrixSuffix -o $outFolder/{}_$matrixPrefix{}${matrixSuffix}_${start_pos}_${end_pos}.png -k $colToSkip -b ${binSizeKb}000 -c {} -s $start_pos -e $end_pos" -- ${all_chromos[@]}
        parallel -i -j $maxJobs -l $maxLoad sh -c "Rscript $script_name  -m $matrixFolder/$matrixPrefix{}$matrixSuffix -o $outFolder/{}_$matrixPrefix{}${matrixSuffix}_${start_pos}_${end_pos}.png -k $colToSkip -b ${binSizeKb}000 -c {} -s $start_pos -e $end_pos" -- ${all_chromos[@]}

    done # end iterating over positions to plot



done # end iterating over datasets


