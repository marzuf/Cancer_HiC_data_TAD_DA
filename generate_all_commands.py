#! /usr/bin/python3



all_data = {
"ENCSR079VIJ_G401": [[],["TCGAkich_norm_kich"]],
"ENCSR401TBQ_Caki2": [[],["TCGAkich_norm_kich"]],
"ENCSR079VIJ_G401ENCSR401TBQ_Caki2": [[],["TCGAkich_norm_kich"]],

"ENCSR312KHQ_SK-MEL-5": [[],["TCGAskcm_lowInf_highInf", "TCGAskcm_wt_mutBRAF", "TCGAskcm_wt_mutCTNNB1"]],
"ENCSR862OGI_RPMI-7951": [[],["TCGAskcm_lowInf_highInf", "TCGAskcm_wt_mutBRAF", "TCGAskcm_wt_mutCTNNB1"]],
"ENCSR312KHQ_SK-MEL-5ENCSR862OGI_RPMI-7951": [[],["TCGAskcm_lowInf_highInf", "TCGAskcm_wt_mutBRAF", "TCGAskcm_wt_mutCTNNB1"]],
  
"ENCSR346DCU_LNCaP": [[],["TCGAprad_norm_prad"]],
"GSE73782_PC3": [[],["TCGAprad_norm_prad"]],
"GSE73782_PC3_ICE": [[],["TCGAprad_norm_prad"]],
"ENCSR346DCU_LNCaPGSE73782_PC3": [[],["TCGAprad_norm_prad"]],
"ENCSR346DCU_LNCaPGSE73782_PC3_ICE": [[],["TCGAprad_norm_prad"]],

"ENCSR444WCZ_A549": [[],["TCGAluad_mutKRAS_mutEGFR", "TCGAluad_nonsmoker_smoker", "TCGAluad_wt_mutKRAS", "TCGAlusc_norm_lusc"]],
"NCI-H460": [[],["TCGAluad_mutKRAS_mutEGFR", "TCGAluad_nonsmoker_smoker", "TCGAluad_wt_mutKRAS", "TCGAlusc_norm_lusc"]],
"ENCSR444WCZ_A549NCI-H460": [[],["TCGAluad_mutKRAS_mutEGFR", "TCGAluad_nonsmoker_smoker", "TCGAluad_wt_mutKRAS", "TCGAlusc_norm_lusc"]],

"ENCSR346DCU_LNCaP": [[],["TCGAprad_norm_prad"]],

"MCF-7":  [[],["TCGAbrca_lum_bas"]],
"ENCSR549MGQ_T47D":  [[],["TCGAbrca_lum_bas"]],
"MCF-7ENCSR549MGQ_T47D":  [[],["TCGAbrca_lum_bas"]],

"Panc1_rep12":  [[],["TCGApaad_wt_mutKRAS"]],

"K562":  [[],["TCGAlaml_wt_mutFLT3"]],

"GSE105318_DLD1":  [[],["TCGAcoad_msi_mss"]],

"GSE105194_spinal_cord":  [[],["TCGAgbm_classical_mesenchymal"]],
"GSE105194_spinal_cord": [[], ["TCGAgbm_classical_proneural"]],
"GSE105194_spinal_cord":  [[],["TCGAgbm_classical_neural"]],
"GSE105194_spinal_cord":  [[],["TCGAlgg_IDHwt_IDHmutnc"]]
}

mystring=""

for hicds in all_data.keys():

	#print(hicds + ":\t" + ','.join(all_data[hicds]) + "\n")
	mystring += hicds + ":\t" + ','.join(all_data[hicds][1]) + "\n"


	hicds_resol = hicds + "_40kb"

	# if this is a consensus
	if len(all_data[hicds][0]) > 0:

		mystring += "Rscript 2_all_chromo_find_consensusTADs.R " + ' '.join(all_data[hicds][0]) + "\n"


	else:



		if hicds == "MCF-7" or hicds == "K562":
		
			mystring += "./aws1_all_chromo_dwd_juicer_hic.sh " + hicds + "\n"

			mystring += "./aws2_all_chromo_hicData_rebinning.sh " + hicds + "\n"

			mystring += "./aws3_all_chromo_prepFile.sh " + hicds + "\n"

			mystring += "./aws4_all_chromo_pre2hic2matrix.sh " + hicds + "\n"

			mystring += "./1_all_chromo_hic_to_TADs.sh " + hicds + "# (set STEP1=0 !!!)	\n"

		
	
		elif hicds == "GSE73782_PC3":
			mystring += "./pc3_all_chromo_dwd_juicer_hic.sh " + hicds + "\n"
			# take the 40kb count files and generate pre file 
			mystring += "./aws3_all_chromo_prepFile.sh " + hicds + "\n"
			# juicer command to generate hic file with norm from pre
			mystring += "./aws4_all_chromo_pre2hic2matrix.sh " + hicds + "\n"
			# hic to TADs for GSE73782_PC3 (set STEP1=0 !!!)
			mystring += "./1_all_chromo_hic_to_TADs.sh " + hicds + "# (set STEP1=0 !!!)	\n"
	
		elif hicds == "GSE73782_PC3_ICE":
			# get count matrix from Yuanlong folder 40kb
			mystring += "./pc3_all_chromo_dwd_juicer_hic.sh " + hicds + "\n"
			# ICE normalize the 40kb counts
			mystring += "./pc3_all_chromo_ice_normalize.sh " + hicds + "\n"
			# hic to TADs for GSE73782_PC3 (set STEP1=0 AND STEP2=0  !!!) ## CHANGE HARD-CODED NORM TO ICE !!!
			mystring += "./1_all_chromo_hic_to_TADs.sh " + hicds + "# (set STEP1=0 and STEP2=0 !!!)	\n"

		else:
			mystring += "./1_all_chromo_hic_to_TADs.sh " + hicds + "\n"

	mystring += "./3_assign_genes.sh " + hicds + "\n"
	mystring += "Rscript check_matResol.R " + hicds_resol + "\n"


	for exprds in all_data[hicds][1]:
	

		mystring += "./run_pipeline.sh " + hicds_resol + " " + exprds + "\n"
		
		mystring += "Rscript check_TCGA_fpkm.R " + hicds_resol + " " + exprds + "\n"

		mystring += "Rscript check_step17_files.R " + hicds_resol + " " + exprds + "\n"

		mystring += "Rscript compare_pipeline_signif.R " + exprds + " " + hicds + "\n"

		mystring += "Rscript compare_pipeline_ranks.R " + exprds + " " + hicds + "\n"

		mystring += "Rscript create_coexpr_sortNoDup_otherTADfile.R " + hicds + " " + exprds + "\n"




	# need to be run after pipeline run (pipeline_geneList.Rdata)

	mystring += "Rscript create_dist_sortNoDup_otherTADfile.R " + hicds_resol + "\n"

	mystring += "Rscript prep_gene_families_TAD_data_otherTADfile.R " + hicds_resol + "\n"
	mystring += "Rscript create_sameFamily_sortNoDup_otherFamFile.R " + hicds_resol + "\n"

	mystring += "Rscript create_sameTAD_sortNoDup_otherTADfile.R " + hicds_resol + "\n" 


print(mystring)