# ./run_v2_pipeline.sh

# retrieve the datasets that have already been run

# ./run_pipeline.sh GSE105381_HepG2_40kb TCGAlihc_norm_lihc 

# TAD_DE_script="./zzz_run_given_step_given_data_v2.sh"

pipScript="./run_pipeline.sh"


# 2v2, 5v2, 6v2, 7v2, 9v2, 11v2

all_hicds=($(ls "PIPELINE/OUTPUT_FOLDER"))

for hicds in "${all_hicds[@]}"; do
	echo $hicds
	all_exprds=($(ls "PIPELINE/OUTPUT_FOLDER/$hicds"))
	
# 	settingFolder="PIPELINE/INPUT_FILES/$hicds"
	
	
	for exprds in "${all_exprds[@]}"; do
	
# 		settingFile="$settingFolder/run_settings_${exprds}.R"

		$pipScript $hicds $exprds	

		#exit 0

		
	done
	
done




