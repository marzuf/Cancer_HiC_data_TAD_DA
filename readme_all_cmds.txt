DATA ORGANIZATION:

PC3_rep_12_100kb
	-> NORM_MAT
	-> TopDom_MAT
	-> TopDom_FILES
	-> FINAL_DOMAINS
	-> genes2tad

#**************************************************************************************
# SCRIPT 1: FROM Hi-C TO TAD CALLING
#**************************************************************************************

./1_hic_to_TADs.sh <cell_line_name> <chromo>
# e.g. ./1_hic_to_TAD.sh PC3_rep12 chr1

# this creates the following files:
PC3_rep_12_100kb/NORM_MAT/PC3_rep_12_chr1_KR_100kb.hic.matrix (step 1)
PC3_rep_12_100kb/TopDom_MAT/PC3_rep_12_chr1_KR_100kb.hic.TopDom.matrix (step 2)
PC3_rep_12_100kb/TopDom_FILES/PC3_rep_12_chr1_KR_100kb.bed (step 3)
PC3_rep_12_100kb/TopDom_FILES/PC3_rep_12_chr1_KR_100kb_final_domains.txt (step 4)

#**************************************************************************************
# SCRIPT 2: CONSENSUS TADs [optional]
#**************************************************************************************

Rscript find_consensusTADs.R  MCF-7 T47D


#**************************************************************************************
# SCRIPT 3: ASSIGN GENES TO NEW TAD LIST
#**************************************************************************************

./3_all_chromo_assign_genes.sh MCF-7


#**************************************************************************************
# PRE-PIPELINE ANALYSIS A: COMPARE THE HI-C DATASETS
#**************************************************************************************


cmp_datasets_matching.R
cmp_datasets_MoC.R
cmp_datasets_nbrTADs.R
cmp_datasets_resol.R


#**************************************************************************************
# SCRIPT 4: RUN THE PIPELINE (prepare setting file + launch the pipeline)
#**************************************************************************************


# ./run_pipeline.sh GSE105194_ENCFF027IEO_astroCerebellum_vs_GSE105957_ENCFF715HDW_astroSpinal TCGAgbm_classical_mesenchymal


#**************************************************************************************
# POST-PIPELINE ANALYSIS A: COMPUTE AUC COEXPR DIST
#**************************************************************************************



all_prep_data_for_AUC_coexprDist_noDist.sh
create_sameTAD_sortNoDup_otherTADfile.R
create_sameFamily_sortNoDup_otherFamFile.R

all_prep_data_for_AUC_coexprDist_onlyDist.sh
create_dist_sortNoDup_otherTADfile.R

all_prep_data_for_AUC_coexprDist_coexpr.sh
create_coexpr_sortNoDup_otherTADfile.R


all_AUC_coexprDist_sortNoDup_otherTADfile_otherFamFile.sh
AUC_coexprDist_withFam_sortNoDup_otherTADfile_otherFamFile.R



#**************************************************************************************
# POST-PIPELINE ANALYSIS B: COMPARE PIPELINE RESULTS CONSENSUS VS. TISSUE-SPECIFIC
#**************************************************************************************


cmp_auc_coexprDist_pipeline_tissue.R

cmp_auc_fcc_coexprDist_pipeline_tissue.R

cmp_auc_fcc_pipeline_tissue.R


cmp_tissue_consensus_fam_TADs_otherFamFile.R


compare_pipeline_results.R
compare_pipeline_signif.R
