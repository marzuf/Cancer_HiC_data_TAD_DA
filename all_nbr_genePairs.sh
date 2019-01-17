#!/usr/bin/bash

# ./all_nbr_genePairs.sh

# Rscript nbr_genePairs_by_dist.R ENCSR079VIJ_G401_40kb TCGAkich_norm_kich

start_time=$(date -R)   

scriptGenePairs="nbr_genePairs_by_dist.R"

# should have hicds + expr ds !!!

all_TAD_files_ds=(
# DS launched 17.01.2018:
"ENCSR079VIJ_G401_40kb TCGAkich_norm_kich"
"ENCSR401TBQ_Caki2_40kb TCGAkich_norm_kich"
"ENCSR079VIJ_G401ENCSR401TBQ_Caki2_40kb TCGAkich_norm_kich"

"ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_lowInf_highInf"
"ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_wt_mutBRAF"
"ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_wt_mutCTNNB1"

"ENCSR312KHQ_SK-MEL-5ENCSR862OGI_RPMI-7951_40kb TCGAskcm_lowInf_highInf"
"ENCSR312KHQ_SK-MEL-5ENCSR862OGI_RPMI-7951_40kb TCGAskcm_wt_mutBRAF"
"ENCSR312KHQ_SK-MEL-5ENCSR862OGI_RPMI-7951_40kb TCGAskcm_wt_mutCTNNB1"

"ENCSR862OGI_RPMI-7951_40kb TCGAskcm_lowInf_highInf"
"ENCSR862OGI_RPMI-7951_40kb TCGAskcm_wt_mutBRAF"
"ENCSR862OGI_RPMI-7951_40kb TCGAskcm_wt_mutCTNNB1"

"ENCSR346DCU_LNCaP_40kb TCGAprad_norm_prad"

"ENCSR444WCZ_A549_40kb TCGAluad_mutKRAS_mutEGFR"
"ENCSR444WCZ_A549_40kb TCGAluad_nonsmoker_smoker"
"ENCSR444WCZ_A549_40kb TCGAluad_wt_mutKRAS"
"ENCSR444WCZ_A549_40kb TCGAlusc_norm_lusc"

"NCI-H460_40kb TCGAluad_mutKRAS_mutEGFR"
"NCI-H460_40kb TCGAluad_nonsmoker_smoker"
"NCI-H460_40kb TCGAluad_wt_mutKRAS"
"NCI-H460_40kb TCGAlusc_norm_lusc"

"ENCSR444WCZ_A549NCI-H460_40kb TCGAluad_mutKRAS_mutEGFR"
"ENCSR444WCZ_A549NCI-H460_40kb TCGAluad_nonsmoker_smoker"
"ENCSR444WCZ_A549NCI-H460_40kb TCGAluad_wt_mutKRAS"
"ENCSR444WCZ_A549NCI-H460_40kb TCGAlusc_norm_lusc"

"K562_40kb  TCGAlaml_wt_mutFLT3"

"MCF-7_40kb TCGAbrca_lum_bas"

"ENCSR549MGQ_T47D_40kb TCGAbrca_lum_bas"

"MCF-7ENCSR549MGQ_T47D_40kb TCGAbrca_lum_bas"

"Panc1_rep12_40kb TCGApaad_wt_mutKRAS"

# missing:
# prostate 2/2; prostate consensus; colorectal; liver; astrocytes 1/2; astrocytes 2/2; astrocyte consensus

)

for ds in "${all_TAD_files_ds[@]}"; do
    
    echo Rscript $scriptGenePairs $ds
    Rscript $scriptGenePairs $ds
done


###################################################################################################################################################
########## END ####################################################################################################################################
###################################################################################################################################################

echo "*** DONE"
end_time=$(date -R)    
echo $start_time
echo $end_time
exit 0

