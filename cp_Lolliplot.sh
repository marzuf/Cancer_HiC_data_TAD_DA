#!/usr/bin/bash

# ./cp_Lolliplot.sh


all_datasets=($( ls PIPELINE/OUTPUT_FOLDER/*/*/*/combined_pval_top20.svg ))

outFold="TOP20_LOLLIPLOTS"
mkdir -p $outFold

for ds in "${all_datasets[@]}"; do

dir1=`dirname $ds`
dir2=`dirname $dir1`

exprds=`basename $dir2`

dir3=`dirname $dir2`
hicds=`basename $dir3`

dsFile=`basename $ds`

#echo $exprds
#echo $hicds

cmd="cp $ds $outFold/${exprds}_${hicds}_${dsFile}"

echo $cmd
$cmd

done


#inFile="PIPELINE/OUTPUT_FOLDER/MCF-7_40kb/TCGAbrca_lum_bas/13_plotTopWilcoxTopCombined/combined_pval_top20.svg"