
######################################################################################################################################################################################################
###################################################################################################################################################################################################### CANCER SUBTYPES (vector)
######################################################################################################################################################################################################

cancer_subAnnot <- c(
  "TCGAstad_msi_gs" = "subtypes",
  "GSE79209_dysp_nodysp" = "lesions",
  "GSE58135_tripleNeg_adjTripleNeg" = "vs_normal",
  "GSE40419_normal_cancer" = "vs_normal",
  "TCGAstad_EBVneg_EBVpos" = "subtypes",
  "GSE77509_normal_ptt" = "vs_normal",
  "GSE77509_normal_tumor" = "vs_normal",
  "GSE77314_normal_tumor" = "vs_normal",
  "GSE71119_dediffSM_MFSM" = "subtypes",
  "GSE58135_ERpos_adjERpos" = "vs_normal",
  "GSE58135_ERpos_tripleNeg" = "subtypes",
  "GSE71119_undiffSM_LMSM" = "subtypes",
  "TCGAlaml_laml_mutFLT3" = "mutation",
  "TCGAthca_thca_mutBRAF" = "mutation",
  "TCGAlihc_lihc_mutCTNNB1" = "mutation",
  "GSE77509_ptt_tumor" = "vs_other",
  "GSE87340_ad_nl" = "vs_normal",
  "GSE81089_normal_nsclc" = "vs_normal",
  "TCGApaad_paad_mutKRAS" = "mutation",
  "TCGAucec_msi_cnl" = "subtypes",        
  "TCGAbrca_lum_bas" = "subtypes",
  "TCGAskcm_skcm_mutBRAF" = "mutation",
  "TCGAluad_luad_mutKRAS" = "mutation",
  "TCGAskcm_skcm_mutCTNNB1" = "mutation",        
  "TCGAacc_acc_mutCTNNB1" = "mutation",
  "GSE74927_neg_pos" = "subtypes",
  "TCGAcrc_msi_mss" =  "subtypes",
  "GSE102073_stic_nostic" = "lesions",
  
  "TCGAcesc_adeno_squam" = "subtypes",
  "TCGAhnsc_HPVneg_HPVpos" = "subtypes",
  "TCGAlgg_IDHwt_IDHmutnc" = "subtypes",
  "TCGAsarc_ddlps_lms" = "subtypes",
  "TCGAsarc_ddlps_mfs" = "subtypes",
  "TCGAsarc_lms_mfs" = "subtypes",
  "TCGAtgct_sem_nonsem" = "subtypes",
  "TCGAcoad_msi_mss" = "subtypes",
  
  "TCGAskcm_lowinf_highInf" = "subtypes",
  "TCGAluad_nonsmoker_smoker" = "subtypes",
  
  "TCGAblca_norm_blca" = "vs_normal",
  "TCGAkich_norm_kich" = "vs_normal",
  "TCGAlusc_norm_lusc" = "vs_normal",
  "TCGAstad_norm_gs" = "vs_normal",
  
  "TCGAprad_norm_prad" = "vs_normal",
  
  "TCGAacc_wt_mutCTNNB1" = "mutation",
  "TCGAgbm_classical_mesenchymal"="subtypes",
  "TCGAgbm_classical_neural"="subtypes",
  "TCGAgbm_classical_proneural"  ="subtypes",
  "TCGAlaml_wt_mutFLT3" = "mutation",
  "TCGAlihc_wt_mutCTNNB1"         = "mutation",
  "TCGAluad_mutKRAS_mutEGFR"="subtypes",
  "TCGAluad_wt_mutKRAS"           = "mutation",
  "TCGApaad_wt_mutKRAS" = "mutation",
  "TCGAskcm_lowInf_highInf"    ="subtypes",  
  "TCGAskcm_wt_mutBRAF" = "mutation",
  "TCGAskcm_wt_mutCTNNB1" = "mutation",        
  "TCGAstad_EBVpos_gs"="subtypes",
  "TCGAthca_mut.RAS_mutBRAF"   =  "subtypes",
  "TCGAthca_wt_mutBRAF"           = "mutation"
  
  
  
)

cancer_subColors <- c(
  subtypes = "green4",
  lesions = "chocolate1",
  vs_normal = "blue3",
  vs_other = "orchid1",
  mutation = "red"
)

stopifnot(cancer_subAnnot %in% names(cancer_subColors))

# to get the color: cancer_subColors[cancer_subAnnot[curr_ds]]