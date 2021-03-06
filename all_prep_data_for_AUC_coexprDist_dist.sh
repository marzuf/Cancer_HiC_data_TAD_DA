#!/usr/bin/bash

# ./all_prep_data_for_AUC_coexprDist_dist.sh

start_time=$(date -R)   

scriptPrepDist="create_dist_sortNoDup_otherTADfile.R"

# DS launched 15.01.2018:
all_TAD_files_ds=(									# => r ### => r 19.01
# "MCF-7_40kb"  # run separately for test
# "ENCSR549MGQ_T47D_40kb"
# "MCF-7ENCSR549MGQ_T47D_40kb"
# "K562_40kb"
# "ENCSR444WCZ_A549_40kb"
# "NCI-H460_40kb"
# "ENCSR444WCZ_A549NCI-H460_40kb"
# "ENCSR346DCU_LNCaP_40kb"
# "Panc1_rep12_40kb"
# "ENCSR079VIJ_G401_40kb"
# "ENCSR401TBQ_Caki2_40kb"
# "ENCSR079VIJ_G401ENCSR401TBQ_Caki2_40kb"
# "ENCSR312KHQ_SK-MEL-5_40kb"
# "ENCSR862OGI_RPMI-7951_40kb"
# "ENCSR312KHQ_SK-MEL-5ENCSR862OGI_RPMI-7951_40kb"
# "GSM2334834_U266_HindIII_40kb"
# "GSM2334832_RPMI-8226_HindIII_40kb"
# "GSM2334834_U266_HindIIIGSM2334832_RPMI-8226_HindIII_40kb"
# "ENCSR834DXR_SK-N-MC_40kb"
# "ENCSR105KFX_SK-N-DZ_40kb"

# added 19.01.2018
# "GSE73782_PC3_40kb"
# "GSE73782_PC3_ICE_40kb"
# "GSE105318_DLD1_40kb"
# "GSE105194_spinal_cord_40kb"
# "ENCSR346DCU_LNCaPGSE73782_PC3_40kb"
# "ENCSR346DCU_LNCaPGSE73782_PC3_ICE_40kb"

# TO RUN 21.01.2018
"GSE105194_cerebellum_40kb"
"GSE105194_spinal_cordGSE105194_cerebellum_40kb"
)
# missing:
# prostate 2/2; prostate consensus; colorectal; liver; astrocytes 1/2; astrocytes 2/2; astrocyte consensus

for ds in "${all_TAD_files_ds[@]}"; do

    echo Rscript $scriptPrepDist $ds
    Rscript $scriptPrepDist $ds
done


###################################################################################################################################################
########## END ####################################################################################################################################
###################################################################################################################################################

echo "*** DONE"
end_time=$(date -R)    
echo $start_time
echo $end_time
exit 0

