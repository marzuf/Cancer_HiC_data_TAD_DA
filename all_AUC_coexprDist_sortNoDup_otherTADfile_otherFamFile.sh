#!/usr/bin/bash

# ./all_AUC_coexprDist_sortNoDup_otherTADfile_otherFamFile.sh

start_time=$(date -R)   

script_name="AUC_coexprDist_withFam_sortNoDup_otherTADfile_otherFamFile.R"

# Rscript AUC_coexprDist_withFam_sortNoDup_otherTADfile_otherFamFile.R ENCSR079VIJ_G401_40kb TCGAkich_norm_kich hgnc
# Rscript AUC_coexprDist_withFam_sortNoDup_otherTADfile_otherFamFile.R ENCSR079VIJ_G401_40kb TCGAkich_norm_kich 
 
#############################################################
all_data=(
# DS launched 16.01.2018:
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

#############################################################
for data in "${all_data[@]}"; do
    echo "> START for $data"
	echo Rscript $script_name $data
	Rscript $script_name $data
done




###################################################################################################################################################
########## END ####################################################################################################################################
###################################################################################################################################################

echo "*** DONE"
end_time=$(date -R)    
echo $start_time
echo $end_time
exit 0

