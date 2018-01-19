#!/usr/bin/env bash

# raw data files location
rawDataDir=/GWD/bioinfo/projects/RD-PTS-EpinovaData-analysis/RNAseq/03trimming
# results files location
resultsDir=/GWD/bioinfo/projects/RD-PTS-EpinovaData-analysis/RNAseq/04qualityAssessment

cd $rawDataDir
# run fastqc
for sample in Sample_RA-*
do
	fastqc --outdir=$resultsDir --casava --threads 8 $sample/RA-*.fastq.gz |& tee -a $resultsDir/fastqc.log.txt
done
