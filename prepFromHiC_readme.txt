!!! BE CAREFUL !!! NEED TO BE DONE INTERCHROMO !!!


0) download the raw data at 10 kb
[command line: juicer dump extract 3 column format 10 kb matrix]
# all_chromo_dwd_juicer_hic.sh #
-> MCF-7/RAW_10kb/MC-7_<chromo>_RAW_10kb.hic.counts

1) 10 kb -> 40 kb raw data using *prepFromHiC_aggregCounts.R*
all_chromo_hicData_rebinning.sh
-> MCF-7/RAW_40kb/MC-7_<chromo>_RAW_40kb.hic.counts

1b) aggregate using HiTC and compare the result with *cmp_rebinning_custom_HiTC.R*




2) prepare a pre.txt file from the rebinned data to create a hic file using *prepFromHic_prepFile.R*
# all_chromo_prepFile.sh #

MCF-7/PRE_40kb/MC-7_chr21_RAW_40kb.pre

3) juicer command to generate hic file from 40kb data
MCF-7/RAW_HIC_INTRA_40kb/MC-7_chr21_intra_40kb.hic  
+ juicer command to extract 40kb KR normalized data
MCF-7_40kb/NORM_MAT/MC-7_chr21_KR_40kb.hic.matrix

all_chromo_pre2hic2matrix.sh


