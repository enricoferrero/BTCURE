#!/usr/bin/env bash

bam_dir="/GWD/bioinfo/projects/RD-PTS-EpinovaData-analysis/RNAseq/06postAlignment"
reference_dir="rsem_reference/GRCh38"
output_dir="rsem_quantification"

cd $bam_dir
for sampleBam in Sample_RA-*; do
	sampleBam_dir=$output_dir/$sampleBam
	mkdir $sampleBam_dir
	rsem-calculate-expression --bam --no-bam-output -p 6 --paired-end --forward-prob 0 $sampleBam $sampleBam_dir/rsem >& $sampleBam_dir/rsem.log

