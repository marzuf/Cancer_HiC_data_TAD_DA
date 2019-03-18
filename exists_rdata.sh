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
	"NCI-H460_40kb"
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
		fi
		
		############################################################################

    done
done


echo "!!! NBR MISSING FILES: $missing"

