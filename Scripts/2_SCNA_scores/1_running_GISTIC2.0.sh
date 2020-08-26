#!/bin/bash

### Script for running GISTIC2.0
# seg files (input files) were downloaded from TCGA using TCGAbiolinks
# then the files were filtered to keep only samples that also presented in ASCAT output and those with age information

cd /opt/GISTIC/

TumourPath=/home/nis/kasit/work/Age_differences_cancer/Data/TCGA_segfiles/Seg_file_Tumour

seg_files=$(find "$TumourPath" -name '*_noCNV_seg_tumour.txt')
marker_file=/mnt/work/larder/kasit/Age_differences_cancer/Data/GISTIC2_hg19_support/genome.info.6.0_hg19.na31_minus_frequent_nan_probes_sorted_2.1.txt
CNV_file=/mnt/work/larder/kasit/Age_differences_cancer/Data/GISTIC2_hg19_support/CNV.hg19.bypos.111213.txt

for f in $seg_files;
do
	base=$(basename $f "_noCNV_seg_tumour.txt")

	echo "Working on: $base"

	# create result directory
	mkdir /home/nis/kasit/work/Age_differences_cancer/Analysis_results/CNAs/1_GISTIC2_CNAs/${base}_GISTIC_mkCNV
	outdir=/home/nis/kasit/work/Age_differences_cancer/Analysis_results/CNAs/1_GISTIC2_CNAs/${base}_GISTIC_mkCNV
	
	# GISTIC2
	/opt/GISTIC/gistic2 \
	-seg $f \
	-cnv $CNV_file \
	-mk $marker_file \
	-refgene /opt/GISTIC/refgenefiles/hg19.UCSC.add_miR.140312.refgene.mat \
	-b $outdir \
	-genegistic 1 \
	-smallmem 1 \
	-qvt 0.25 \
	-ta 0.25 \
	-td 0.25 \
	-broad 1 \
	-brlen 0.7 \
	-conf 0.95 \
	-armpeel 1 \
	-savegene 1

	echo "Finish: $base"
done;

