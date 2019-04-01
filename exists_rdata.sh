#!/usr/bin/bash
# ./exists_rdata.sh 10b_runEmpPvalProdSignedRatio/emp_pval_prodSignedRatio.Rdata

rdatafile="$1"

all_hic_dataset=(
	"ENCSR444WCZ_A549_40kb"
	"NCI-H460_40kb"
	"ENCSR444WCZ_A549NCI-H460_40kb"

	"ENCSR079VIJ_G401_40kb"
	"ENCSR401TBQ_Caki2_40kb"
	"ENCSR079VIJ_G401ENCSR401TBQ_Caki2_40kb"

	"ENCSR312KHQ_SK-MEL-5_40kb"
	"ENCSR862OGI_RPMI-7951_40kb"
	"ENCSR312KHQ_SK-MEL-5ENCSR862OGI_RPMI-7951_40kb"

	"ENCSR549MGQ_T47D_40kb"
	"MCF-7_40kb"
	"MCF-7ENCSR549MGQ_T47D_40kb"

	"ENCSR346DCU_LNCaP_40kb"
	"GSE73782_PC3_40kb"
	"GSE73782_PC3_ICE_40kb"
	"ENCSR346DCU_LNCaPGSE73782_PC3_40kb"
	"ENCSR346DCU_LNCaPGSE73782_PC3_ICE_40kb"

	"GSE105194_cerebellum_40kb"
	"GSE105194_spinal_cord_40kb"
	"GSE105194_spinal_cordGSE105194_cerebellum_40kb"

	"GSE105381_HepG2_40kb"

 	"ENCSR346DCU_LNCaP_40kb"

	"GSE105318_DLD1_40kb"

	"Panc1_rep12_40kb"

	"K562_40kb" #-> 6v2 running
)

i=0
missing=0

for hic_dataset in "${all_hic_dataset[@]}"; do

	#echo "> START hic_dataset = $hic_dataset"

	if [[ "$hic_dataset" == "ENCSR444WCZ_A549_40kb" || "$hic_dataset" == "NCI-H460_40kb" || "$hic_dataset" == "ENCSR444WCZ_A549NCI-H460_40kb" ]]; then
		all_expr_dataset=(
		"TCGAluad_norm_luad"
		"TCGAluad_mutKRAS_mutEGFR"
		"TCGAluad_nonsmoker_smoker"
		"TCGAluad_wt_mutKRAS"
		"TCGAlusc_norm_lusc"
		)
	elif [[ "$hic_dataset" == "MCF-7_40kb" || "$hic_dataset" == "ENCSR549MGQ_T47D_40kb" || "$hic_dataset" == "MCF-7ENCSR549MGQ_T47D_40kb" ]]; then
		all_expr_dataset=(
			"TCGAbrca_lum_bas"
		)
	elif [[ "$hic_dataset" == "GSE105381_HepG2_40kb" ]]; then
		all_expr_dataset=(
			"TCGAlihc_norm_lihc"
			"TCGAlihc_wt_mutCTNNB1"
		)
	elif [[ "$hic_dataset" == "K562_40kb" ]]; then
		all_expr_dataset=(
			"TCGAlaml_wt_mutFLT3"
		)
	elif [[ "$hic_dataset" == "GSE105318_DLD1_40kb" ]]; then
		all_expr_dataset=(
			"TCGAcoad_msi_mss"
		)
	elif [[ "$hic_dataset" == "Panc1_rep12_40kb" ]]; then
		all_expr_dataset=(
			"TCGApaad_wt_mutKRAS"
		)
	elif [[ "$hic_dataset" == "GSE73782_PC3_40kb" || "$hic_dataset" == "GSE73782_PC3_ICE_40kb" || "$hic_dataset" == "ENCSR346DCU_LNCaP_40kb" || "$hic_dataset" == "ENCSR346DCU_LNCaPGSE73782_PC3_40kb" || "$hic_dataset" == "ENCSR346DCU_LNCaPGSE73782_PC3_ICE_40kb" ]]; then

		all_expr_dataset=(
			"TCGAprad_norm_prad"
		)

	elif [[ "$hic_dataset" == "ENCSR079VIJ_G401_40kb" || "$hic_dataset" == "ENCSR401TBQ_Caki2_40kb" || "$hic_dataset" == "ENCSR079VIJ_G401ENCSR401TBQ_Caki2_40kb" ]]; then

		all_expr_dataset=(
			"TCGAkich_norm_kich"
		)
	elif [[ "$hic_dataset" == "ENCSR312KHQ_SK-MEL-5_40kb" || "$hic_dataset" == "ENCSR862OGI_RPMI-7951_40kb" || "$hic_dataset" == "ENCSR312KHQ_SK-MEL-5ENCSR862OGI_RPMI-7951_40kb" ]]; then

		all_expr_dataset=(
			"TCGAskcm_lowInf_highInf"
			"TCGAskcm_wt_mutBRAF"
			"TCGAskcm_wt_mutCTNNB1"
		)
	elif [[ "$hic_dataset" == "Panc1_rep12_40kb" ]]; then
		all_expr_dataset=(
				"TCGAprad_norm_prad"
		)
	elif [[ "$hic_dataset" == "GSE105194_cerebellum_40kb" || "$hic_dataset" == "GSE105194_spinal_cord_40kb" || "$hic_dataset" == "GSE105194_spinal_cordGSE105194_cerebellum_40kb" ]]; then
		all_expr_dataset=(
			"TCGAgbm_classical_mesenchymal"
			"TCGAgbm_classical_proneural"
			"TCGAgbm_classical_neural"
			"TCGAlgg_IDHwt_IDHmutnc"
		)
	elif [[ "$hic_dataset" == "GSE105194_cerebellum_40kb" || "$hic_dataset" == "GSE105194_spinal_cord_40kb" || "$hic_dataset" == "GSE105194_spinal_cordGSE105194_cerebellum_40kb" ]]; then
		all_expr_dataset=(
			"TCGAgbm_classical_mesenchymal"
			"TCGAgbm_classical_proneural"
			"TCGAgbm_classical_neural"
			"TCGAlgg_IDHwt_IDHmutnc"
		)
	else
	    echo "invalid hicds"
	    exit 1
	fi

	for expr_dataset in "${all_expr_dataset[@]}"; do

		((i++))

		#echo "> START hic_dataset = $hic_dataset - expr_dataset = $expr_dataset (i=$i)"


		checkDir="PIPELINE/OUTPUT_FOLDER/$hic_dataset/$expr_dataset"

		if [[ ! -d $checkDir ]]; then
		    echo "invalid $checkDir"
		    exit 1
		fi
		
		filename="$checkDir/$rdatafile"

		if [[ ! -f  $filename ]]; then
		    echo "> START hic_dataset = $hic_dataset - expr_dataset = $expr_dataset (i=$i)"
		    echo "... $filename does not exist !"
		    #exit 1
		    ((missing++))
        #else 
            #echo "-- > $filename - OK! "
		fi
		
		############################################################################

    done
done


echo "!!! NBR MISSING FILES: $missing"

############################## # 23.03.2019
#./exists_rdata.sh 0_prepGeneData/pipeline_geneList.Rdata                           # 0
#./exists_rdata.sh 7v2_runPermutationsMeanTADCorr/meanCorr_permDT.Rdata             # 11
#./exists_rdata.sh 10b_runEmpPvalProdSignedRatio/emp_pval_prodSignedRatio.Rdata     # 0
#./exists_rdata.sh 1_runGeneDE/DE_topTable.Rdata                                    # 0
#./exists_rdata.sh 10_runEmpPvalMeanTADCorr/emp_pval_meanCorr.Rdata                 # 0
#./exists_rdata.sh 10v2_runEmpPvalMeanTADCorr/emp_pval_meanCorr.Rdata               # 0 (30.03)                      # 53 (23.03)         
#./exists_rdata.sh 11_runEmpPvalCombined/emp_pval_combined.Rdata                    # 0
#./exists_rdata.sh 170revision2EZH2_score_auc_pval_permGenes/auc_ratios.Rdata       # 0
#./exists_rdata.sh 2_runWilcoxonTAD/wilcox_ttest_meanTAD_qq.Rdata                   # 0
#./exists_rdata.sh 2v2_runWilcoxonTAD/wilcox_pairedTAD_meanExpr_wilcoxStat.Rdata    # 0
#./exists_rdata.sh 3_runMeanTADLogFC/all_meanLogFC_TAD.Rdata                        # 0
#./exists_rdata.sh   4_runMeanTADCorr/all_meanCorr_TAD.Rdata                        # 0
#./exists_rdata.sh 4v2_runConcordTADCorr/all_concordCorr_TAD.Rdata                  # 0
#./exists_rdata.sh 4v3_runConcordTADCorr/all_concordCorr_TAD.Rdata                  # 0
#./exists_rdata.sh 4vAll_runConcordTADCorr/all_concordCorr_TAD.Rdata                # 0
#./exists_rdata.sh 5_runPermutationsMedian/permutationsDT.Rdata                     # 0
#./exists_rdata.sh 6_runPermutationsMeanLogFC/meanLogFC_permDT.Rdata                # 0
#./exists_rdata.sh 6v2onlyW_runPermutationsWilcoxStat/wilcoxStat_permDT.Rdata       # 0
#./exists_rdata.sh 6v2_runPermutationsWilcoxStat/*DT.Rdata                          # 55 -> ok
#./exists_rdata.sh 7_runPermutationsMeanTADCorr/meanCorr_permDT.Rdata               # 0
#./exists_rdata.sh 7v2_runPermutationsMeanTADCorr/meanCorr_permDT.Rdata             # 0 (30.03)        # 11 (23.03) -> missing 5v2 and 7v2, running
#./exists_rdata.sh 8c_runAllDown/prodSignedRatio_permDT.Rdata                       # 0
#./exists_rdata.sh 9_runEmpPvalMeanTADLogFC/emp_pval_meanLogFC.Rdata                # 0

#./exists_rdata.sh 5v2_runPermutations/permutationsDT.Rdata
#./exists_rdata.sh 6v2onlyW_runPermutationsWilcoxStat/wilcoxStat_permDT.Rdata      
#./exists_rdata.sh 7v2_runPermutationsMeanTADCorr/meanCorr_permDT.Rdata            
#./exists_rdata.sh 9v2_runEmpPvalWilcoxStat/emp_pval_wilcoxStat.Rdata
#./exists_rdata.sh 10v2_runEmpPvalMeanTADCorr/emp_pval_meanCorr.Rdata              
#./exists_rdata.sh 10b_runEmpPvalProdSignedRatio/emp_pval_prodSignedRatio.Rdata    


############################################################################
############################################################################
############################################################################
#./exists_rdata.sh 7v2_runPermutationsMeanTADCorr/meanCorr_permDT.Rdata 
#> START hic_dataset = GSE105194_spinal_cord_40kb - expr_dataset = TCGAgbm_classical_neural (i=44)
#... PIPELINE/OUTPUT_FOLDER/GSE105194_spinal_cord_40kb/TCGAgbm_classical_neural/7v2_runPermutationsMeanTADCorr/meanCorr_permDT.Rdata does not exist !
#> START hic_dataset = GSE105194_spinal_cord_40kb - expr_dataset = TCGAlgg_IDHwt_IDHmutnc (i=45)
#... PIPELINE/OUTPUT_FOLDER/GSE105194_spinal_cord_40kb/TCGAlgg_IDHwt_IDHmutnc/7v2_runPermutationsMeanTADCorr/meanCorr_permDT.Rdata does not exist !
#> START hic_dataset = GSE105194_spinal_cordGSE105194_cerebellum_40kb - expr_dataset = TCGAgbm_classical_mesenchymal (i=46)
#... PIPELINE/OUTPUT_FOLDER/GSE105194_spinal_cordGSE105194_cerebellum_40kb/TCGAgbm_classical_mesenchymal/7v2_runPermutationsMeanTADCorr/meanCorr_permDT.Rdata does not exist !
#> START hic_dataset = GSE105194_spinal_cordGSE105194_cerebellum_40kb - expr_dataset = TCGAgbm_classical_proneural (i=47)
#... PIPELINE/OUTPUT_FOLDER/GSE105194_spinal_cordGSE105194_cerebellum_40kb/TCGAgbm_classical_proneural/7v2_runPermutationsMeanTADCorr/meanCorr_permDT.Rdata does not exist !
#> START hic_dataset = GSE105194_spinal_cordGSE105194_cerebellum_40kb - expr_dataset = TCGAgbm_classical_neural (i=48)
#... PIPELINE/OUTPUT_FOLDER/GSE105194_spinal_cordGSE105194_cerebellum_40kb/TCGAgbm_classical_neural/7v2_runPermutationsMeanTADCorr/meanCorr_permDT.Rdata does not exist !
#> START hic_dataset = GSE105194_spinal_cordGSE105194_cerebellum_40kb - expr_dataset = TCGAlgg_IDHwt_IDHmutnc (i=49)
#... PIPELINE/OUTPUT_FOLDER/GSE105194_spinal_cordGSE105194_cerebellum_40kb/TCGAlgg_IDHwt_IDHmutnc/7v2_runPermutationsMeanTADCorr/meanCorr_permDT.Rdata does not exist !
#> START hic_dataset = GSE105381_HepG2_40kb - expr_dataset = TCGAlihc_norm_lihc (i=50)
#... PIPELINE/OUTPUT_FOLDER/GSE105381_HepG2_40kb/TCGAlihc_norm_lihc/7v2_runPermutationsMeanTADCorr/meanCorr_permDT.Rdata does not exist !
#> START hic_dataset = GSE105381_HepG2_40kb - expr_dataset = TCGAlihc_wt_mutCTNNB1 (i=51)
#... PIPELINE/OUTPUT_FOLDER/GSE105381_HepG2_40kb/TCGAlihc_wt_mutCTNNB1/7v2_runPermutationsMeanTADCorr/meanCorr_permDT.Rdata does not exist !
#> START hic_dataset = GSE105318_DLD1_40kb - expr_dataset = TCGAcoad_msi_mss (i=53)
#... PIPELINE/OUTPUT_FOLDER/GSE105318_DLD1_40kb/TCGAcoad_msi_mss/7v2_runPermutationsMeanTADCorr/meanCorr_permDT.Rdata does not exist !
#> START hic_dataset = Panc1_rep12_40kb - expr_dataset = TCGApaad_wt_mutKRAS (i=54)
#... PIPELINE/OUTPUT_FOLDER/Panc1_rep12_40kb/TCGApaad_wt_mutKRAS/7v2_runPermutationsMeanTADCorr/meanCorr_permDT.Rdata does not exist !
#> START hic_dataset = K562_40kb - expr_dataset = TCGAlaml_wt_mutFLT3 (i=55)
#... PIPELINE/OUTPUT_FOLDER/K562_40kb/TCGAlaml_wt_mutFLT3/7v2_runPermutationsMeanTADCorr/meanCorr_permDT.Rdata does not exist !
#!!! NBR MISSING FILES: 11

#PIPELINE/OUTPUT_FOLDER/ENCSR079VIJ_G401_40kb/TCGAkich_norm_kich/10v2_runEmpPvalMeanTADCorr/emp_pval_meanCorr.Rdata
#PIPELINE/OUTPUT_FOLDER/ENCSR444WCZ_A549_40kb/TCGAluad_mutKRAS_mutEGFR/10v2_runEmpPvalMeanTADCorr/emp_pval_meanCorr.Rdata

