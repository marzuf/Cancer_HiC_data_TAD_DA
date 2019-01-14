#!/usr/bin/bash

# ./all_prep_data_for_AUC_coexprDist_noDist.sh

start_time=$(date -R)   

scriptSameTAD="create_sameTAD_sortNoDup_otherTADfile.R"


all_TAD_files_ds=(
"GSE105194_ENCFF027IEO_astroCerebellum_vs_GSE105957_ENCFF715HDW_astroSpinal"
"GSE105318_ENCFF439QFU_DLD1"
"GSE105465_ENCFF777DUA_Caki2_vs_GSE105235_ENCFF235TGH_G401"
"GSE105566_ENCFF358MNA_Panc1"
"GSE105600_ENCFF852YOE_A549_vs_GSE105725_ENCFF697NNX_NCIH460"
"GSE106022_ENCFF614EKT_RPMI7951_vs_GSE105491_ENCFF458OWO_SKMEL5"
"GSE63525_K562"
"GSM1631185_MCF7_vs_GSE75070_MCF7_shGFP"
"GSM1631185_MCF7_vs_GSE75070_MCF7_shGFP_vs_GSE105697_ENCFF364CWZ_T47D"
)


for ds in "${all_TAD_files_ds[@]}"; do

    echo Rscript $scriptSameTAD $ds
    Rscript $scriptSameTAD $ds

done


###################################################################################################################################################
########## END ####################################################################################################################################
###################################################################################################################################################

echo "*** DONE"
end_time=$(date -R)    
echo $start_time
echo $end_time
exit 0

