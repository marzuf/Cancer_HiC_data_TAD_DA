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

maxJobs=40
maxLoad=70

script_name="draw_matrix.R"

binSizeKb=40
colToSkip=3

#all_chromos=( "chr1" "chr9" "chr21" )
#all_chromos=( "chr5" )
all_chromos=( "chr"{1..22} ) 

all_dataset_id=(
## breast
"GSE105697_ENCFF364CWZ_T47D"
"GSM1631185_MCF7"
"GSM1631185_MCF7"
## lung
"GSE105600_ENCFF852YOE_A549"
"GSE105725_ENCFF697NNX_NCIH460"
## kidney
"GSE105465_ENCFF777DUA_Caki2"
"GSE105235_ENCFF235TGH_G401"
# leukemia
"GSE63525_K562"
# prostate
"GSE105557_ENCFF270HJX_LNCaP"
# pancreas
"GSE105566_ENCFF358MNA_Panc1"
# colorectal
"GSE105318_ENCFF439QFU_DLD1"
"GSE105318_ENCFF714TMN_DLD1_int"
"GSE105318_ENCFF714TMN_DLD1_intICE"
# skin
"GSE106022_ENCFF614EKT_RPMI7951"
"GSE105491_ENCFF458OWO_SKMEL5"
# astrocyte cerebellum
"GSE105194_ENCFF027IEO_ASTROCEREB"
"GSE105194_ENCFF122YID_ASTROCEREB_int"
"GSE105194_ENCFF122YID_ASTROCEREB_intICE"
# astrocyte spinal cord
"GSE105957_ENCFF478UBU_ASTROSPINAL"
"GSE105957_ENCFF715HDW_ASTROSPINAL_int"
"GSE105957_ENCFF715HDW_ASTROSPINAL_intICE"
)


for dataID in "${all_dataset_id[@]}"; do

    ################
    ######## BREAST
    ################

    if [[ $dataID == "GSE105697_ENCFF364CWZ_T47D" ]]; then
        matrixFolder="/mnt/etemp/marie/Dixon2018_integrative_data/breast/T47D/ENCSR549MGQ_GSE105697"
        matrixPrefix="GSE105697_ENCFF364CWZ_chromatin_interactions_hg19_"
        matrixSuffix="_TopDom.matrix"
        outFolder="DRAW_ENCODE_MATRIX/$dataID"
    elif [[ $dataID == "GSM1631185_MCF7" ]]; then
        matrixFolder="/mnt/etemp/marie/Dixon2018_integrative_data/breast/MCF7/GSM1631185_GSE66733/Hi-C_MCF7_MCF10A_processed_HiCfiles/Heatmaps/chrxchr/40kb/MCF7"
        matrixPrefix="HiCStein-MCF7-WT__hg19__"
        matrixSuffix="__C-40000-iced_TopDom.matrix"
        outFolder="DRAW_ENCODE_MATRIX/GSM1631185_MCF7"
    elif [[ $dataID == "GSM1631185_MCF7" ]]; then
        matrixFolder="/mnt/etemp/marie/Dixon2018_integrative_data/breast/MCF7/GSM1942100_GSM1942101_GSE75070/shGFP"
        matrixPrefix="GSE75070_HiCStein-MCF7-shGFP_hg19_"
        matrixSuffix="_C-40000-iced_TopDom.matrix"
        outFolder="DRAW_ENCODE_MATRIX/GSE75070_MCF7_shGFP"

    ################
    ### LUNG
    ################
    elif [[ $dataID == "GSE105600_ENCFF852YOE_A549" ]]; then
        matrixFolder="/mnt/etemp/marie/Dixon2018_integrative_data/lung/A549/ENCSR444WCZ"
        matrixPrefix="GSE105600_ENCFF852YOE_chromatin_interactions_hg19_"
        matrixSuffix="_TopDom.matrix"
        outFolder="DRAW_ENCODE_MATRIX/GSE105600_ENCFF852YOE_A549"

    elif [[ $dataID == "GSE105725_ENCFF697NNX_NCIH460" ]]; then
        matrixFolder="/mnt/etemp/marie/Dixon2018_integrative_data/lung/NCI-H460/ENCSR489OCU"
        matrixPrefix="GSE105725_ENCFF697NNX_chromatin_interactions_hg19_"
        matrixSuffix="_TopDom.matrix"
        outFolder="DRAW_ENCODE_MATRIX/GSE105725_ENCFF697NNX_NCIH460"

    ################
    ### PANCREAS
    ################
    elif [[ $dataID == "GSE105566_ENCFF358MNA_Panc1" ]]; then
        matrixFolder="/mnt/etemp/marie/Dixon2018_integrative_data/pancreas/Panc1/GSE105566/ENCFF358MNA"
        matrixPrefix="GSE105566_ENCFF358MNA_chromatin_interactions_hg19_"
        matrixSuffix="_TopDom.matrix"
        outFolder="DRAW_ENCODE_MATRIX/GSE105566_ENCFF358MNA_Panc1"

    ################
    ### KIDNEY
    ################
    elif [[ $dataID == "GSE105465_ENCFF777DUA_Caki2" ]]; then
        matrixFolder="/mnt/etemp/marie/Dixon2018_integrative_data/kidney/Caki2/GSE105465/ENCFF777DUA"
        matrixPrefix="GSE105465_ENCFF777DUA_chromatin_interactions_hg19_"
        matrixSuffix="_TopDom.matrix"
        outFolder="DRAW_ENCODE_MATRIX/GSE105465_ENCFF777DUA_Caki2"

    elif [[ $dataID == "GSE105235_ENCFF235TGH_G401" ]]; then
        matrixFolder="/mnt/etemp/marie/Dixon2018_integrative_data/kidney/G401/GSE105235/ENCFF235TGH"
        matrixPrefix="GSE105235_ENCFF235TGH_chromatin_interactions_hg19_"
        matrixSuffix="_TopDom.matrix"
        outFolder="DRAW_ENCODE_MATRIX/GSE105235_ENCFF235TGH_G401"

    ################
    ### COLORECTAL
    ################
    elif [[ $dataID == "GSE105318_ENCFF439QFU_DLD1" ]]; then
        matrixFolder="/mnt/etemp/marie/Dixon2018_integrative_data/colon/DLD1/GSE105318/ENCFF439QFU"
        matrixPrefix="GSE105318_ENCFF439QFU_chromatin_interactions_hg19_"
        matrixSuffix="_TopDom.matrix"
        outFolder="DRAW_ENCODE_MATRIX/GSE105318_ENCFF439QFU_DLD1"

    elif [[ $dataID == "GSE105318_ENCFF714TMN_DLD1_int" ]]; then
        matrixFolder="/mnt/etemp/marie/Dixon2018_integrative_data/colon/DLD1/GSE105318/ENCFF714TMN"
        matrixPrefix="GSE105318_ENCFF714TMN_chromatin_interactions_hg19_"
        matrixSuffix="_TopDom.matrix"
        outFolder="DRAW_ENCODE_MATRIX/GSE105318_ENCFF714TMN_DLD1_int"

    elif [[ $dataID == "GSE105318_ENCFF714TMN_DLD1_intICE" ]]; then
        matrixFolder="/mnt/etemp/marie/Dixon2018_integrative_data/colon/DLD1/GSE105318/ENCFF714TMN"
        matrixPrefix="GSE105318_ENCFF714TMN_chromatin_interactions_hg19_"
        matrixSuffix="_ICE_TopDom.matrix"
        outFolder="DRAW_ENCODE_MATRIX/GSE105318_ENCFF714TMN_DLD1_intICE"

    ################
    ### SKIN
    ################
    elif [[ $dataID == "GSE106022_ENCFF614EKT_RPMI7951" ]]; then
        matrixFolder="/mnt/etemp/marie/Dixon2018_integrative_data/skin/RPMI-7951/GSE106022/ENCFF614EKT"
        matrixPrefix="GSE106022_ENCFF614EKT_chromatin_interactions_hg19_"
        matrixSuffix="_TopDom.matrix"
        outFolder="DRAW_ENCODE_MATRIX/GSE106022_ENCFF614EKT_RPMI7951"

    elif [[ $dataID == "GSE105491_ENCFF458OWO_SKMEL5" ]]; then
        matrixFolder="/mnt/etemp/marie/Dixon2018_integrative_data/skin/SK-MEL-5/GSE105491/ENCFF458OWO"
        matrixPrefix="GSE105491_ENCFF458OWO_chromatin_interactions_hg19_"
        matrixSuffix="_TopDom.matrix"
        outFolder="DRAW_ENCODE_MATRIX/GSE105491_ENCFF458OWO_SKMEL5"

    ################
    ### LYMPHOCYTE
    ################
    elif [[ $dataID == "GSE63525_K562" ]]; then
        matrixFolder="/mnt/etemp/marie/Dixon2018_integrative_data/leukemia/K562/GSE63525"
        matrixPrefix="GSE63525_K562_40kb_ICE_"
        matrixSuffix="_TopDom.matrix"
        outFolder="DRAW_ENCODE_MATRIX/GSE63525_K562"

    ################
    ### ASTROCYTE - CEREBELLUM
    ################
    elif [[ $dataID == "GSE105194_ENCFF027IEO_ASTROCEREB" ]]; then
        matrixFolder="/mnt/etemp/marie/Dixon2018_integrative_data/astrocyte/cerebellum/GSE105194/ENCFF027IEO"
        matrixPrefix="GSE105194_ENCFF027IEO_chromatin_interactions_hg19_"
        matrixSuffix="_TopDom.matrix"
        outFolder="DRAW_ENCODE_MATRIX/GSE105194_ENCFF027IEO_ASTROCEREB"

    elif [[ $dataID == "GSE105194_ENCFF122YID_ASTROCEREB_int" ]]; then
        matrixFolder="/mnt/etemp/marie/Dixon2018_integrative_data/astrocyte/cerebellum/GSE105194/ENCFF122YID"
        matrixPrefix="GSE105194_ENCFF122YID_chromatin_interactions_hg19_"
        matrixSuffix="_TopDom.matrix"
        outFolder="DRAW_ENCODE_MATRIX/GSE105194_ENCFF122YID_ASTROCEREB_int"

    elif [[ $dataID == "GSE105194_ENCFF122YID_ASTROCEREB_intICE" ]]; then
        matrixFolder="/mnt/etemp/marie/Dixon2018_integrative_data/astrocyte/cerebellum/GSE105194/ENCFF122YID"
        matrixPrefix="GSE105194_ENCFF122YID_chromatin_interactions_hg19_"
        matrixSuffix="_ICE_TopDom.matrix"
        outFolder="DRAW_ENCODE_MATRIX/GSE105194_ENCFF122YID_ASTROCEREB_intICE"

    ################
    ### ASTROCYTE - SPINAL CORD
    ################
    elif [[ $dataID == "GSE105957_ENCFF715HDW_ASTROSPINAL_int" ]]; then
        matrixFolder="/mnt/etemp/marie/Dixon2018_integrative_data/astrocyte/spinal_cord/GSE105957/ENCFF715HDW"
        matrixPrefix="GSE105957_ENCFF715HDW_chromatin_interactions_hg19_"
        matrixSuffix="_TopDom.matrix"
        outFolder="DRAW_ENCODE_MATRIX/GSE105957_ENCFF715HDW_ASTROSPINAL"

    elif [[ $dataID == "GSE105957_ENCFF478UBU_ASTROSPINAL" ]]; then
        matrixFolder="/mnt/etemp/marie/Dixon2018_integrative_data/astrocyte/spinal_cord/GSE105957/ENCFF478UBU"
        matrixPrefix="GSE105957_ENCFF478UBU_chromatin_interactions_hg19_"
        matrixSuffix="_TopDom.matrix"
        outFolder="DRAW_ENCODE_MATRIX/GSE105957_ENCFF478UBU_ASTROSPINAL_int"

    elif [[ $dataID == "GSE105957_ENCFF478UBU_ASTROSPINAL_intICE" ]]; then
        matrixFolder="/mnt/etemp/marie/Dixon2018_integrative_data/astrocyte/spinal_cord/GSE105957/ENCFF478UBU"
        matrixPrefix="GSE105957_ENCFF715HDW_chromatin_interactions_hg19_"
        matrixSuffix="_ICE_TopDom.matrix"
        outFolder="DRAW_ENCODE_MATRIX/GSE105957_ENCFF478UBU_ASTROSPINAL_intICE"


    ################
    ### PROSTATE
    ################
    elif [[ $dataID == "GSE105557_ENCFF270HJX_LNCaP" ]]; then
        matrixFolder="/mnt/etemp/marie/Dixon2018_integrative_data/prostate/LNCaP/GSE105557/ENCFF270HJX"
        matrixPrefix="GSE105557_ENCFF270HJX_chromatin_interactions_hg19_"
        matrixSuffix="_TopDom.matrix"
        outFolder="DRAW_ENCODE_MATRIX/GSE105557_ENCFF270HJX_LNCaP"

    else
        echo "invalid dataID"
        exit 0
    fi

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


