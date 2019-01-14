
topdomPattern <- "_final_domains.txt$"

cl_labs <- c(
"MCF-7" = "JR",
"ENCSR549MGQ_T47D" = "JD",

"DLD1" = "JD",

"ENCSR079VIJ_G401" = "JD",
"ENCSR401TBQ_Caki2" = "JD",

"HepG2" = "JD",

"ENCSR444WCZ_A549" = "JD",
"NCI-H460" = "JD",

"K562" = "ELA",

"GSM2334834_U266_HindIII" = "CL",
"GSM2334832_RPMI-8226_HindIII" ="CL",

"ENCSR834DXR_SK-N-MC" = "JD",
"ENCSR105KFX_SK-N-DZ" = "JD",

"Panc1_rep12" = "JD",

"ENCSR346DCU_LNCaP" = "JD",
"GSE73782_PC3" = "JAK",

"ENCSR312KHQ_SK-MEL-5" = "JD",
"ENCSR862OGI_RPMI-7951" = "JD",

"pipelineConsensus" = NA
)
stopifnot(!duplicated(names(cl_labs)))


cl_names <- c(
"breastCl1" = "MCF-7",
"breastCl2" = "ENCSR549MGQ_T47D",

"colorectalCl1" = "DLD1",                     # !!! RUNNING !!!

"kidneyCl1" = "ENCSR079VIJ_G401",
"kidneyCl2" = "ENCSR401TBQ_Caki2",

"liver" = "HepG2",							# !!! RUNNING !!!

"lungCl1" = "ENCSR444WCZ_A549",
"lungCl2" = "NCI-H460",

"lymphoblast" = "K562",

"myelomaCl1" = "GSM2334834_U266_HindIII",
"myelomaCl2" = "GSM2334832_RPMI-8226_HindIII",

"neuroblastomaCl1" = "ENCSR834DXR_SK-N-MC",
"neuroblastomaCl2" = "ENCSR105KFX_SK-N-DZ",

"pancreasCl1" = "Panc1_rep12",

"prostateCl1" = "ENCSR346DCU_LNCaP",
"prostateCl2" = "GSE73782_PC3",				# !!! RUNNING !!!

"skinCl1" = "ENCSR312KHQ_SK-MEL-5",
"skinCl2" = "ENCSR862OGI_RPMI-7951",

"pipelineConsensus" = "pipelineConsensus"
)

stopifnot( length(cl_names) == length(cl_labs))
stopifnot( cl_names %in% names(cl_labs))


# do not take neuroblastoma consensus -> derived from different metastatic sites

cl_names["breastConsensus"] <- paste0(cl_names["breastCl1"], cl_names["breastCl2"]) 
cl_names["kidneyConsensus"] <- paste0(cl_names["kidneyCl1"], cl_names["kidneyCl2"]) 
cl_names["lungConsensus"] <- paste0(cl_names["lungCl1"], cl_names["lungCl2"]) 
cl_names["myelomaConsensus"] <- paste0(cl_names["myelomaCl1"], cl_names["myelomaCl2"]) 
cl_names["prostateConsensus"] <- paste0(cl_names["prostateCl1"], cl_names["prostateCl2"]) 		# !!! RUNNING !!!
cl_names["skinConsensus"] <- paste0(cl_names["skinCl1"], cl_names["skinCl2"]) 
	
stopifnot(!duplicated(names(cl_names)))
stopifnot(!duplicated(cl_names))


