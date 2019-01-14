
topdomPattern <- "_final_domains.txt$"




cl_names <- c(
"breastCl1" = "MCF-7",
"breastCl2" = "ENCSR549MGQ_T47D",

"lungCl1" = "A549",
"lungCl2" = "NCI-H460",

"lymphoblast" = "K562",

"pancreasCl1" = "Panc1",

"prostateCl1" = "PC3",
"prostateCl2" = "LNCaP",

"kidneyCl1" = "Caki2",
"kidneyCl2" = "G401",

"skinCl1" = "RPMI-7951",
"skinCl2" = "SK-MEL-5"


)

cl_names["breastConsensus"] <- paste0(cl_names["breastCl1"], cl_names["breastCl2"]) 
cl_names["kidneyConsensus"] <- paste0(cl_names["kidneyCl1"], cl_names["kidneyCl2"]) 
cl_names["lungConsensus"] <- paste0(cl_names["lungCl1"], cl_names["lungCl2"]) 
cl_names["prostateConsensus"] <- paste0(cl_names["prostateCl1"], cl_names["prostateCl2"]) 
cl_names["skinConsensus"] <- paste0(cl_names["skinCl1"], cl_names["skinCl2"]) 


cl_labs <- c(
"MCF-7" = "barutcu",
"ENCSR549MGQ_T47D" = "?",

"A549" = "?",
"NCI-H460" = "?",

"K562" = "rao",

"Panc1" = "?",

"PC3" = "?",
"LNCaP" = "?",

"Caki2" ="?",
"G401" ="?",

"RPMI-7951" = "?",
"SK-MEL-5" = "?"


)


#================================================
#================================================ BREAST
#================================================

breastConsensusname <-  "GSM1631185_MCF7_vs_GSE75070_MCF7_shGFP_vs_GSE105697_ENCFF364CWZ_T47D"
breastConsensusFold <- file.path("FIND_CONSENSUS_TADS", breastConsensusname)
breastConsensusFiles <- list.files(breastConsensusFold, full.names=T, pattern = consensusPattern)
stopifnot(length(breastConsensusFiles) > 0)
breast_consensus_chromos <- unique(gsub("(chr.+)_conservedTADs.txt", "\\1", basename(breastConsensusFiles)))
stopifnot(length(breast_consensus_chromos) > 0)

mcf7Consensusname <-  "GSM1631185_MCF7_vs_GSE75070_MCF7_shGFP"
mcf7ConsensusFold <- file.path("FIND_CONSENSUS_TADS",mcf7Consensusname)
mcf7ConsensusFiles <- list.files(mcf7ConsensusFold, full.names=T, pattern = consensusPattern)
stopifnot(length(mcf7ConsensusFiles) > 0)
mcf7_consensus_chromos <- unique(gsub("(chr.+)_conservedTADs.txt", "\\1", basename(mcf7ConsensusFiles)))
stopifnot(length(mcf7_consensus_chromos) > 0)

# breast/MCF7/GSM1631185_GSE66733/MCF7/TopDom/GSM1631185_MCF7_40kb_chr1_final_domains.txt
breastCL1name <- "GSM1631185_GSE66733"
breastFold1 <- file.path("breast/MCF7/", breastCL1name, "MCF7/TopDom")
breastCL1Files <- list.files(breastFold1, full.names=T, pattern = topdomPattern)
stopifnot(length(breastCL1Files) > 0)
breast1_chromos <- unique(gsub(".+(chr.+)_final_domains.txt$", "\\1", basename(breastCL1Files)))
stopifnot(length(breast1_chromos) > 0)

# breast/MCF7/GSM1942100_GSM1942101_GSE75070/shGFP/TopDom/GSE75070_MCF7_shGFP_40kb_chr1_final_domains.txt
breastCL2name <- "GSM1942100_GSM1942101_GSE75070"
breastFold2 <- file.path("breast/MCF7", breastCL2name, "shGFP/TopDom")
breastCL2Files <- list.files(breastFold2, full.names=T, pattern = topdomPattern)
stopifnot(length(breastCL2Files) > 0)
breast2_chromos <- unique(gsub(".+(chr.+)_final_domains.txt$", "\\1", basename(breastCL2Files)))
stopifnot(length(breast2_chromos) > 0)

# breast/T47D/ENCSR549MGQ_GSE105697/TopDom/GSE105697_ENCFF364CWZ_T47D_40kb_chr1_final_domains.txt
breastCL3name <- "ENCSR549MGQ_GSE105697"
breastFold3 <- file.path("breast/T47D", breastCL3name, "TopDom")
breastCL3Files <- list.files(breastFold3, full.names=T, pattern = topdomPattern)
stopifnot(length(breastCL3Files) > 0)
breast3_chromos <- unique(gsub(".+(chr.+)_final_domains.txt$", "\\1", basename(breastCL3Files)))
stopifnot(length(breast3_chromos) > 0)

#================================================
#================================================ LUNG
#================================================

lungConsensusname <-  "GSE105600_ENCFF852YOE_A549_vs_GSE105725_ENCFF697NNX_NCIH460"
lungConsensusFold <- file.path("FIND_CONSENSUS_TADS", lungConsensusname)
lungConsensusFiles <- list.files(lungConsensusFold, full.names=T, pattern = consensusPattern)
stopifnot(length(lungConsensusFiles) > 0)
lung_consensus_chromos <- unique(gsub("(chr.+)_conservedTADs.txt", "\\1", basename(lungConsensusFiles)))
stopifnot(length(lung_consensus_chromos) > 0)

# lung/A549/ENCSR444WCZ/TopDom/GSE105600_ENCFF852YOE_A549_40kb_chr1_final_domains.txt
lungCL1name <- "ENCSR444WCZ"
lungFold1 <- file.path("lung/A549", lungCL1name, "TopDom")
lungCL1Files <- list.files(lungFold1, full.names=T, pattern = topdomPattern)
stopifnot(length(lungCL1Files) > 0)
lung1_chromos <- unique(gsub(".+(chr.+)_final_domains.txt$", "\\1", basename(lungCL1Files)))
stopifnot(length(lung1_chromos) > 0)

# lung/NCI-H460/ENCSR489OCU/TopDom/GSE105725_ENCFF697NNX_NCIH460_40kb_chr1_final_domains.txt
lungCL2name <- "ENCSR489OCU"
lungFold2 <- file.path("lung/NCI-H460", lungCL2name, "TopDom")
lungCL2Files <- list.files(lungFold2, full.names=T, pattern = topdomPattern)
stopifnot(length(lungCL2Files) > 0)
lung2_chromos <- unique(gsub(".+(chr.+)_final_domains.txt$", "\\1", basename(lungCL2Files)))
stopifnot(length(lung2_chromos) > 0)

#================================================
#================================================ PANCREAS
#================================================


#pancreas/Panc1/GSE105566/ENCFF358MNA/TopDom GSE105566_ENCFF358MNA_Panc1_40kb_chr.+final_domains.txt 
pancreasCL1name <- "ENCFF358MNA"
pancreasFold1 <- file.path("pancreas/Panc1/GSE105566", pancreasCL1name, "TopDom")
pancreasCL1Files <- list.files(pancreasFold1, full.names=T, pattern = topdomPattern)
stopifnot(length(pancreasCL1Files) > 0)
pancreas1_chromos <- unique(gsub(".+(chr.+)_final_domains.txt$", "\\1", basename(pancreasCL1Files)))
stopifnot(length(pancreas1_chromos) > 0)


#prostate/LNCaP/GSE105557/ENCFF270HJX/TopDom GSE105557_ENCFF270HJX_LNCaP_40kb_chr.+final_domains.txt 
prostateCL1name <- "ENCFF270HJX"
prostateFold1 <- file.path("prostate/LNCaP/GSE105557", prostateCL1name, "TopDom")
prostateCL1Files <- list.files(prostateFold1, full.names=T, pattern = topdomPattern)
stopifnot(length(prostateCL1Files) > 0)
prostate1_chromos <- unique(gsub(".+(chr.+)_final_domains.txt$", "\\1", basename(prostateCL1Files)))
stopifnot(length(prostate1_chromos) > 0)

#================================================
#================================================ KIDNEY
#================================================

#kidney/Caki2/GSE105465/ENCFF777DUA/TopDom GSE105465_ENCFF777DUA_Caki2_40kb_chr.+final_domains.txt 
kidneyCL1name <- "ENCFF777DUA"
kidneyFold1 <- file.path("kidney/Caki2/GSE105465", kidneyCL1name, "TopDom")
kidneyCL1Files <- list.files(kidneyFold1, full.names=T, pattern = topdomPattern)
stopifnot(length(kidneyCL1Files) > 0)
kidney1_chromos <- unique(gsub(".+(chr.+)_final_domains.txt$", "\\1", basename(kidneyCL1Files)))
stopifnot(length(kidney1_chromos) > 0)

#kidney/G401/GSE105235/ENCFF235TGH/TopDom GSE105235_ENCFF235TGH_G401_40kb_chr.+final_domains.txt 
kidneyCL2name <- "ENCFF235TGH"
kidneyFold2 <- file.path("kidney/G401/GSE105235", kidneyCL2name, "TopDom")
kidneyCL2Files <- list.files(kidneyFold2, full.names=T, pattern = topdomPattern)
stopifnot(length(kidneyCL2Files) > 0)
kidney2_chromos <- unique(gsub(".+(chr.+)_final_domains.txt$", "\\1", basename(kidneyCL2Files)))
stopifnot(length(kidney2_chromos) > 0)

kidneyConsensusname <-  "GSE105465_ENCFF777DUA_Caki2_vs_GSE105235_ENCFF235TGH_G401"
kidneyConsensusFold <- file.path("FIND_CONSENSUS_TADS", kidneyConsensusname)
kidneyConsensusFiles <- list.files(kidneyConsensusFold, full.names=T, pattern = consensusPattern)
stopifnot(length(kidneyConsensusFiles) > 0)
kidney_consensus_chromos <- unique(gsub("(chr.+)_conservedTADs.txt", "\\1", basename(kidneyConsensusFiles)))
stopifnot(length(kidney_consensus_chromos) > 0)


#================================================
#================================================ SKIN
#================================================


#skin/RPMI-7951/GSE106022/ENCFF614EKT/TopDom GSE106022_ENCFF614EKT_RPMI7951_40kb_chr.+final_domains.txt 
skinCL1name <- "ENCFF614EKT"
skinFold1 <- file.path("skin/RPMI-7951/GSE106022", skinCL1name, "TopDom")
skinCL1Files <- list.files(skinFold1, full.names=T, pattern = topdomPattern)
stopifnot(length(skinCL1Files) > 0)
skin1_chromos <- unique(gsub(".+(chr.+)_final_domains.txt$", "\\1", basename(skinCL1Files)))
stopifnot(length(skin1_chromos) > 0)

#skin/SK-MEL-5/GSE105491/ENCFF458OWO/TopDom GSE105491_ENCFF458OWO_SKMEL5_40kb_chr.+final_domains.txt 
skinCL2name <- "ENCFF458OWO"
skinFold2 <- file.path("skin/SK-MEL-5/GSE105491", skinCL2name, "TopDom")
skinCL2Files <- list.files(skinFold2, full.names=T, pattern = topdomPattern)
stopifnot(length(skinCL2Files) > 0)
skin2_chromos <- unique(gsub(".+(chr.+)_final_domains.txt$", "\\1", basename(skinCL2Files)))
stopifnot(length(skin2_chromos) > 0)

skinConsensusname <-  "GSE106022_ENCFF614EKT_RPMI7951_vs_GSE105491_ENCFF458OWO_SKMEL5"
skinConsensusFold <- file.path("FIND_CONSENSUS_TADS", skinConsensusname)
skinConsensusFiles <- list.files(skinConsensusFold, full.names=T, pattern = consensusPattern)
stopifnot(length(skinConsensusFiles) > 0)
skin_consensus_chromos <- unique(gsub("(chr.+)_conservedTADs.txt", "\\1", basename(skinConsensusFiles)))
stopifnot(length(skin_consensus_chromos) > 0)




#================================================
#================================================ COLORECTAL
#================================================


#/mnt/etemp/marie/Dixon2018_integrative_data/colon/DLD1/GSE105318/ENCFF439QFU/GSE105318_ENCFF439QFU_chromatin_interactions_hg19_chr1_TopDom.matrix → float
#/mnt/etemp/marie/Dixon2018_integrative_data/colon/DLD1/GSE105318/ENCFF714TMN/GSE105318_ENCFF714TMN_chromatin_interactions_hg19_chr1_TopDom.matrix → int ???
#"coloCL1" = "GSE105318_ENCFF439QFU",
#"coloCL2" = "GSE105318_ENCFF714TMN", # int



coloCL1name <- "ENCFF439QFU"
coloFold1 <- file.path("colon/DLD1/GSE105318", coloCL1name, "TopDom")
coloCL1Files <- list.files(coloFold1, full.names=T, pattern = topdomPattern)
stopifnot(length(coloCL1Files) > 0)
colo1_chromos <- unique(gsub(".+(chr.+)_final_domains.txt$", "\\1", basename(coloCL1Files)))
stopifnot(length(colo1_chromos) > 0)


## => int matrix
#coloCL2name <- "ENCFF714TMN"
#coloFold2 <- file.path("colon/DLD1/GSE105318", coloCL2name, "TopDom")
#coloCL2Files <- list.files(coloFold2, full.names=T, pattern = topdomPattern)
#stopifnot(length(coloCL2Files) > 0)
#colo2_chromos <- unique(gsub(".+(chr.+)_final_domains.txt$", "\\1", basename(coloCL2Files)))
#stopifnot(length(colo2_chromos) > 0)


                    # NOT DONE
                    #coloConsensusname <-  "GSE106022_ENCFF614EKT_RPMI7951_vs_GSE105491_ENCFF458OWO_SKMEL5"
                    #coloConsensusFold <- file.path("FIND_CONSENSUS_TADS", coloConsensusname)
                    #coloConsensusFiles <- list.files(coloConsensusFold, full.names=T, pattern = consensusPattern)
                    #stopifnot(length(coloConsensusFiles) > 0)
                    #colo_consensus_chromos <- unique(gsub("(chr.+)_conservedTADs.txt", "\\1", basename(coloConsensusFiles)))
                    #stopifnot(length(colo_consensus_chromos) > 0)

#================================================
#================================================ ASTROCYTES
#================================================


#astrocyte/cerebellum/GSE105194/ENCFF027IEO/GSE105194_ENCFF027IEO_chromatin_interactions_hg19_chr1_TopDom.matrix → float
#astrocyte/cerebellum/GSE105194/ENCFF122YID/GSE105194_ENCFF122YID_chromatin_interactions_hg19_chr1_TopDom.matrix → int ???
#astrocyte/spinal_cord/GSE105957/ENCFF478UBU/GSE105957_ENCFF478UBU_chromatin_interactions_hg19_chr1_TopDom.matrix → int ???
#astrocyte/spinal_cord/GSE105957/ENCFF715HDW/GSE105957_ENCFF715HDW_chromatin_interactions_hg19_chr1_TopDom.matrix → float

# "astroCL1" = "GSE105194_ENCFF027IEO",
# "astroCL2"= "GSE105194_ENCFF122YID", # int
# "astroCL3"= "GSE105957_ENCFF715HDW",
#"astroCL4"= "GSE105957_ENCFF478UBU", # int

astroCL1name <- "ENCFF027IEO"
astroFold1 <- file.path("astrocyte/cerebellum/GSE105194", astroCL1name, "TopDom")
astroCL1Files <- list.files(astroFold1, full.names=T, pattern = topdomPattern)
stopifnot(length(astroCL1Files) > 0)
astro1_chromos <- unique(gsub(".+(chr.+)_final_domains.txt$", "\\1", basename(astroCL1Files)))
stopifnot(length(astro1_chromos) > 0)


# int
#astroCL2name <- "ENCFF122YID"
#astroFold2 <- file.path("astrocyte/cerebellum/GSE105194/ENCFF027IEO", astroCL2name, "TopDom")
#astroCL2Files <- list.files(astroFold2, full.names=T, pattern = topdomPattern)
#stopifnot(length(astroCL2Files) > 0)
#astro1_chromos <- unique(gsub(".+(chr.+)_final_domains.txt$", "\\1", basename(astroCL2Files)))
#stopifnot(length(astro2_chromos) > 0)


astroCL3name <- "ENCFF715HDW"
astroFold3 <- file.path("astrocyte/spinal_cord/GSE105957", astroCL3name, "TopDom")
astroCL3Files <- list.files(astroFold3, full.names=T, pattern = topdomPattern)
stopifnot(length(astroCL3Files) > 0)
astro3_chromos <- unique(gsub(".+(chr.+)_final_domains.txt$", "\\1", basename(astroCL3Files)))
stopifnot(length(astro3_chromos) > 0)

# int
#astroCL4name <- "ENCFF478UBU"
#astroFold4 <- file.path("astrocyte/spinal_cord/GSE105957", astroCL4name, "TopDom")
#astroCL4Files <- list.files(astroFold4, full.names=T, pattern = topdomPattern)
#stopifnot(length(astroCL4Files) > 0)
#astro1_chromos <- unique(gsub(".+(chr.+)_final_domains.txt$", "\\1", basename(astroCL4Files)))
#stopifnot(length(astro4_chromos) > 0)


# based only on astrocL1 and astroCL3 (float matrices) FIND_CONSENSUS_TADS/GSE105194_ENCFF027IEO_astroCerebellum_vs_GSE105957_ENCFF715HDW_astroSpinal/chrX_conservedTADs.txt
# NOT DONE
astroConsensusname <-  "GSE105194_ENCFF027IEO_astroCerebellum_vs_GSE105957_ENCFF715HDW_astroSpinal"
astroConsensusFold <- file.path("FIND_CONSENSUS_TADS", astroConsensusname)
astroConsensusFiles <- list.files(astroConsensusFold, full.names=T, pattern = consensusPattern)
stopifnot(length(astroConsensusFiles) > 0)
astro_consensus_chromos <- unique(gsub("(chr.+)_conservedTADs.txt", "\\1", basename(astroConsensusFiles)))
stopifnot(length(astro_consensus_chromos) > 0)


#================================================
#================================================ LIVER
#================================================



#================================================
#================================================ LYMPHOBlAST k562 - LEUKEMIAS
#================================================

#/mnt/etemp/marie/Dixon2018_integrative_data/leukemia/K562/GSE63525/GSE63525_K562_40kb_ICE_chr1_TopDom.matrix
#"lympho1" = "GSE63525_K562"

#lymphoCL1name <- "GSE63525_K562"
#lymphoFold1 <- file.path("leukemia/K562/GSE63525", lymphoCL1name, "TopDom")
#lymphoCL1Files <- list.files(lymphoFold1, full.names=T, pattern = topdomPattern)
#stopifnot(length(lymphoCL1Files) > 0)
#lympho1_chromos <- unique(gsub(".+(chr.+)_final_domains.txt$", "\\1", basename(lymphoCL1Files)))
#stopifnot(length(lympho1_chromos) > 0)


#================================================
#================================================ PIPELINE
#================================================


pipConsensusname <- "pipeline_TopDom"
pipConsensusFold <- file.path(setDir, "/mnt/ed4/marie/TAD_call_pipeline_TopDom", "consensus_TopDom_covThresh_r0.6_t80000_v0_w-1_final")
pipConsensusFiles <- list.files(pipConsensusFold, full.names=T, pattern = consensusPattern)
pip_consensus_chromos <- unique(gsub("(chr.+)_conservedTADs.txt", "\\1", basename(pipConsensusFiles)))


#================================================
#================================================ INTERSECT CHROMOS
#================================================

intersectChromos <- Reduce(intersect, list(
  pip_consensus_chromos,
  breast_consensus_chromos, mcf7_consensus_chromos, 
  lung_consensus_chromos,
  breast1_chromos, breast2_chromos, breast3_chromos,
  lung1_chromos, lung2_chromos,
  pancreas1_chromos,
  prostate1_chromos,
  kidney1_chromos, kidney2_chromos,
  skin1_chromos, skin2_chromos,

#lympho1_chromo,

#colo_consensus_chromos,
colo1_chromos,
#colo2_chromos,

astro_consensus_chromos, # based only on astrocL1 and astroCL3 (float matrices)
astro1_chromos,
#astro2_chromos,
astro3_chromos
#astro4_chromos,



))

#================================================
#================================================ ALL DATASETS VECTOR
#================================================


# 6 consensus + 11 datasets
all_ds <- c(
"pipConsensus",
"breastConsensus", 
"mcf7Consensus", 
"lungConsensus", 
"kidneyConsensus", 
"skinConsensus",


"breastCL1", 
"breastCL2", 
"breastCL3",
"lungCL1", 
"lungCL2",
"pancreasCL1", 
"prostateCL1",
"kidneyCL1", 
"kidneyCL2",
"skinCL1", 
"skinCL2",

#"coloConsensus",
"coloCL1",
#"coloCL2", # int

"astroConsensus",  # based only on astrocL1 and astroCL3 (float matrices)
"astroCL1",
#"astroCL2", # int
"astroCL3"
#"astroCL4",# int

#"lympho1"

)

#================================================
#================================================ NAME SETTINGS
#================================================


ds_mapping <- c(
"breastCL1" = "HiCStein-MCF7-WT__hg19__.+TopDom.matrix" , 
"breastCL2" = "GSE75070_HiCStein-MCF7-shGFP_hg19_.+TopDom.matrix",
"breastCL3" = "GSE105697_ENCFF364CWZ", 

"lungCL1" = "GSE105600_ENCFF852YOE",
"lungCL2" = "GSE105725_ENCFF697NNX",

"pancreasCL1" = "GSE105566_ENCFF358MNA", 

"prostateCL1" = "GSE105557_ENCFF270HJX",

"kidneyCL1" = "GSE105465_ENCFF777DUA", 
"kidneyCL2" =  "GSE105235_ENCFF235TGH",

"skinCL1" = "GSE106022_ENCFF614EKT", 
"skinCL2" = "GSE105491_ENCFF458OWO",


"coloCL1" = "GSE105318_ENCFF439QFU",
"coloCL2" = "GSE105318_ENCFF714TMN", # int

"astroCL1" = "GSE105194_ENCFF027IEO",
 "astroCL2"= "GSE105194_ENCFF122YID", # int
"astroCL3"= "GSE105957_ENCFF715HDW",
"astroCL4"= "GSE105957_ENCFF478UBU", # int

"lympho1" = "GSE63525_K562_40kb_ICE"


)

