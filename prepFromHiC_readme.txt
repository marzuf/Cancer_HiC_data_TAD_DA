!!! BE CAREFUL !!! NEED TO BE DONE INTERCHROMO !!!


0) download the raw data at 10 kb
[command line: juicer dump extract 3 column format 10 kb matrix]

1) 10 kb -> 40 kb raw data using prepFromHiC_aggregCounts.R

1b) aggregate using HiTC and compare the result with cmp_rebinning_custom_HiTC.R

2) prepare a pre.txt file from the rebinned data to create a hic file using prepFromHic_prepNewHic.R
I need the hic file for the 40 kb raw data then for the normalization using juicer

3) KR normalization of the hic that contains 40 kb raw data 
[command line: juicer addNorm]


  addNorm <input_HiC_file> [-w genome-wide-resolution] -F

Required:

    <input_HiC_file>: File to normalize; this will delete any previous normalizations Optional:
    -w <genome-wide resolution>: Smallest resolution to calculate genome-wide resolution; e.g., if 10000, genome-wide normalizations will be calculated for 2.5Mb, 1Mb, 500Kb, 250Kb, 100Kb, 50Kb, 25Kb, and 10Kb but not for 5Kb. Note that genome-wide resolution can be very expensive in terms of memory; this flags allows for a memory/normalization trade-off. If not set or set to 0, no genome-wide resolutions will be calculated
    -F: Do not calculate normalizations for fragment-delimited matrices
    -d: For genome-wide normalization, include intra-chromosomal matrices; by default, inter-only matrices are used.
