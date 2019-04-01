#!/usr/bin/bash

#./run_pipeline_allDS.sh 4v2
#./run_pipeline_allDS.sh 4vAll
#./run_pipeline_allDS.sh 6v2onlyW
#./run_pipeline_allDS.sh 10b
#./run_pipeline_allDS.sh 9v2
#./run_pipeline_allDS.sh 10v2
#./run_pipeline_allDS.sh 7v2
# 5v2: all_hic_dataset running

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

for hic_dataset in "${all_hic_dataset[@]}"; do

	echo "> START hic_dataset = $hic_dataset"

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

		echo "> START hic_dataset = $hic_dataset - expr_dataset = $expr_dataset (i=$i)"


		checkDir="PIPELINE/OUTPUT_FOLDER/$hic_dataset/$expr_dataset"

		if [[ ! -d $checkDir ]]; then
		    echo "invalid $checkDir"
		    exit 1
		fi

			
		start_time=$(date -R)    
		#set -e

		if [[ $# != 1 ]]; then
		    echo "invalid # of arguments"
		    exit 1
		fi
		TAD_DE_pipSteps="$1"

		echo "*** START ***"
		echo "... > Hi-C dataset: $hic_dataset"
		echo "... > Gene expression dataset: $expr_dataset"

		#********************** HARD-CODED SETTINGS FOR THE PIPELINE ********************************************
		step1=1     # prepare setting file
		step2=1    # run the pipeline

		runDir="/mnt/etemp/marie/Cancer_HiC_data_TAD_DA"

		TAD_DE_pipDir="/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom"
		TAD_DE_script="./zzz_run_given_step_given_data_v2.sh"
		old_inputFolder="/mnt/ed4/marie/scripts/TAD_DE_pipeline/SETTING_FILES_cleanInput"
		nPermut="10000"
		#nPermut="10"
		ncpu="50"

		Rexec=`which Rscript`

		new_inputFolder="$runDir/PIPELINE/INPUT_FILES/$hic_dataset"
		mkdir -p $new_inputFolder
		outputFolder="$runDir/PIPELINE/OUTPUT_FOLDER/$hic_dataset/$expr_dataset"
		mkdir -p $outputFolder

		echo "!!! IMPORTANT HARD-CODED SETTINGS !!!"
		echo "... ! step1 = $step1"
		echo "... ! step2 = $step2"

		echo "... ! old_inputFolder = $old_inputFolder"
		echo "... ! nPermut = $nPermut"
		echo "... ! ncpu = $ncpu"

		if [[ $step2 -eq 1 ]]; then
			echo "... ! TAD_DE_pipDir = $TAD_DE_pipDir"
			echo "... ! TAD_DE_script = $TAD_DE_script"
			echo "... ! TAD_DE_pipStep(s): ${TAD_DE_pipSteps[*]}"
		fi

		###################################### FUNCTION DEFINITIONS

		function mvBack {
		  echo "... go back to my folder"
		  cd $runDir  
		}
		trap mvBack EXIT

		runCMD() {
		  echo "> $1"
		  eval $1
		}

		checkFile() {
		if [[ ! -f  $2 ]]; then
		    echo "... $1 ($2) does not exist !"
		    exit 1
		fi
		}
		############################################################################

		old_setting_file="$old_inputFolder/run_settings_${expr_dataset}.R"

		checkFile old_setting_file $old_setting_file

		new_setting_file="$new_inputFolder/run_settings_${expr_dataset}.R"

		runCMD "cp $old_setting_file $new_inputFolder"

		# !!! TO CHECK FORMAT ZZZZZ
		# $hic_dataset/genes2tad/all_assigned_regions.txt
		# MCF-7/genes2tad/all_assigned_regions.txt
		# TADposDT.txt
		#chr1    chr1_TAD1       750001  1300000
		new_TADpos_file="$hic_dataset/genes2tad/all_assigned_regions.txt"

		# !!! TO CHECK FORMAT ZZZZZ
		# $hic_dataset/genes2tad/all_genes_positions.txt
		# MCF-7/genes2tad/all_genes_positions.txt
		#gene2tadDT.txt
		#12893       chr1    761586  762902  chr1_TAD1
		new_gene2tad_file="$hic_dataset/genes2tad/all_genes_positions.txt"

		checkFile new_TADpos_file $new_TADpos_file

		checkFile new_gene2tad_file $new_gene2tad_file

		if [[ "$step1" -eq 1 ]] ; then

			cat >> ${new_setting_file} <<- EOM

					# > file edited: `date -R` 

					# path to output folder:
					pipOutFold <- "${outputFolder}"

					# OVERWRITE THE DEFAULT SETTINGS FOR INPUT FILES - use TADs from the current Hi-C dataset 
					TADpos_file <- paste0(setDir, "`realpath $new_TADpos_file`")
									#chr1    chr1_TAD1       750001  1300000
									#chr1    chr1_TAD2       2750001 3650000
									#chr1    chr1_TAD3       3650001 4150000

					gene2tadDT_file <- paste0(setDir, "`realpath $new_gene2tad_file`")
									#LINC00115       chr1    761586  762902  chr1_TAD1
									#FAM41C  chr1    803451  812283  chr1_TAD1
									#SAMD11  chr1    860260  879955  chr1_TAD1
									#NOC2L   chr1    879584  894689  chr1_TAD1

					# overwrite main_settings.R: nCpu <- 25
					nCpu <- ${ncpu}

					# *************************************************************************************************************************
					# ************************************ SETTINGS FOR PERMUTATIONS (5#_, 8c_)
					# *************************************************************************************************************************

					# number of permutations
					nRandomPermut <- $nPermut
					gene2tadAssignMethod <- "maxOverlap"
					nRandomPermutShuffle <- $nPermut
					step8_for_permutGenes <- TRUE
					step8_for_randomTADsFix <- FALSE
					step8_for_randomTADsGaussian <- FALSE
					step8_for_randomTADsShuffle <- FALSE
					step14_for_randomTADsShuffle <- FALSE

			EOM
			echo "WRITTEN: ${new_setting_file}"


		echo "> END STEP1:" $(date -R)

		checkFile new_setting_file $new_setting_file

		fi # end if STEP1

		###################################################################################################################################################
		### STEP2: MOVE TO THE TAD DE PIPELINE DIRECTORY TO LAUNCH THE TAD DE PIPELINE
		###################################################################################################################################################
		if [[ "$step2" -eq 1 ]] ; then
			echo "> START STEP2:" $(date -R)
			runCMD "cd $TAD_DE_pipDir"

			echo $TAD_DE_script ${new_setting_file} `echo ${TAD_DE_pipSteps[*]}`
			$TAD_DE_script ${new_setting_file} `echo ${TAD_DE_pipSteps[*]}`

			cd $runDir
			echo "> END STEP2:" $(date -R)
		fi




	done #iterate exprds
done #iterate hicds

###################################################################################################################################################
########## END ####################################################################################################################################
###################################################################################################################################################
echo "*** DONE"
end_time=$(date -R)    
echo $start_time
echo $end_time
exit 0

