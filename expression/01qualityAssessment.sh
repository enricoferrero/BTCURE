#!/usr/bin/env bash

# raw data files location
rawDataDir=/GWD/bioinfo/projects/RD-PTS-EpinovaData-raw/RNAseq
cd "$rawDataDir"
# results files location
resultsDir=/GWD/bioinfo/projects/RD-PTS-EpinovaData-analysis/RNAseq/01qualityAssessment

# check md5sum
md5sum -c md5sum.txt | tee "$resultsDir"/md5sum.log.txt
# run fastqc
for sample in Sample_RA-*
do
	fastqc --outdir="$resultsDir" --casava --threads 8 $sample/*.fastq.gz | tee -a "$resultsDir"/fastqc.stdout.txt
done
